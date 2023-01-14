# Adapted from this tutorial: 

# install.packages(c("rayshader", "raster", "sp"))
# if rayshader does not load, download Xquartz and then retry installing rgl
# if raster does not load, must install rgdal, first install Homebrew 
# (https://stackoverflow.com/a/63236256) then do `brew install gdal`

library(rayshader)
library(sp)
library(raster)
library(scales)
library(terra)

# Load aerial image 
# assumes you are in the main directory already 
fn = "stratmap20-nc-cir-12in-caparea_3097253d1_20200107.jp2"
home_r = raster::raster(fn, band = 1)
home_g = raster::raster(fn, band = 2)
home_b = raster::raster(fn, band = 3)

home_rbg = raster::stack(home_r, home_g, home_b)
#raster::plotRGB(home_rbg)

# Load the elevation DEM as a mosaic (not sure which tile I need, load all)
# From TNRIS - https://data.tnris.org/collection?c=447db89a-58ee-4a1b-a61f-b918af2fb0bb&geo=-97.96926882045099,30.522691292708643,-97.96567323076479,30.525645751161917#11.48/30.5334/-97.9671
# https://gis.stackexchange.com/a/408215
efiles <- list.files("stratmap21-28cm-50cm-bexar-travis_3097253_dem", full.names = T, pattern = "*.tif")
ic <- terra::sprc(lapply(efiles, terra::rast))
elevation = terra::mosaic(ic)
elevation2 = raster::raster(elevation)

# Check to make sure the projections match
raster::crs(home_r)
raster::crs(elevation2)

# Project elevation raster to match the aerial image
elevation_utm = raster::projectRaster(elevation2, crs = crs(home_r), method = "bilinear")
raster::crs(elevation_utm) 

# Go to google maps to get these values from the area you want
# Most zoomed in
top_right = c(y=, x=)
bottom_left = c(y=-, x=)

extent_latlong = sp::SpatialPoints(rbind(bottom_left, top_right), proj4string=sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
extent_utm = sp::spTransform(extent_latlong, raster::crs(elevation_utm))

e = raster::extent(extent_utm)

# Now crop it 
home_rgb_cropped = raster::crop(home_rbg, e)
elevation_cropped = raster::crop(elevation_utm, e)

names(home_rgb_cropped) = c("r","g","b")

home_r_cropped = rayshader::raster_to_matrix(home_rgb_cropped$r)
home_g_cropped = rayshader::raster_to_matrix(home_rgb_cropped$g)
home_b_cropped = rayshader::raster_to_matrix(home_rgb_cropped$b)

homeel_matrix = rayshader::raster_to_matrix(elevation_cropped)

home_rgb_array = array(0,dim=c(nrow(home_r_cropped),ncol(home_r_cropped),3))

home_rgb_array[,,1] = home_r_cropped/255 #Red layer
home_rgb_array[,,2] = home_g_cropped/255 #Blue layer
home_rgb_array[,,3] = home_b_cropped/255 #Green layer

home_rgb_array = aperm(home_rgb_array, c(2,1,3))

#plot_map(home_rgb_array)

home_rgb_contrast = scales::rescale(home_rgb_array,to=c(0,1))

plot_map(home_rgb_contrast)

# Plot the 3D Map
plot_3d(home_rgb_contrast, homeel_matrix, windowsize = c(1400,900),
        zoom=0.55, phi=30, theta=115, fov=120)

#rayshader::save_3dprint('home.stl')

render_snapshot(title_text = "Home, Leander, Texas | Imagery: CAPCOG Aerial | DEM: 1ft LIDAR", 
                title_bar_color = "#1f5214", title_color = "white", title_bar_alpha = 1)
