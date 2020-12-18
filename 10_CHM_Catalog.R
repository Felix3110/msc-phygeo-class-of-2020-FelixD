#------------------------------------------------------------------------------
# Type: control script 
# Name: 10_CHM_Catalog.R
# Author: Chris Reudenbach, creuden@gmail.com
# Description:  script creates a canopy height model from generic Lidar 
#              las data using the lidR package. In addition the script cuts the area to 
#              a defined extent using the lidR catalog concept.
#              Furthermore the data is tiled into a handy format even for poor memory 
#              and a canopy height model is calculated
# Data: regular las LiDAR data sets 
# Output: Canopy heightmodel RDS file of the resulting catalog
# Copyright: Chris Reudenbach, Thomas Nauss 2017,2020, GPL (>= 3)
#------------------------------------------------------------------------------

## clean your environment
rm(list=ls()) 

# 0 - load packages
#-----------------------------

library("future")

## dealing with the crs warnings is cumbersome and complex
## you may reduce the warnings with uncommenting  the following line
## for a deeper  however rarely less confusing understanding have a look at:
## https://rgdal.r-forge.r-project.org/articles/CRS_projections_transformations.html
## https://www.r-spatial.org/r/2020/03/17/wkt.html
rgdal::set_thin_PROJ6_warnings(TRUE)


# 1 - source files
#-----------------
source(file.path(envimaR::alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd",
                                          alt_env_id = "COMPUTERNAME",
                                          alt_env_value = "PCRZP",
                                          alt_env_root_folder = "F:/BEN/edu/mpg-envinsys-plygrnd"),
                 "/msc-phygeo-class-of-2020-FelixD/src/000_setup.R"))


# 2 - define variables
#---------------------

## define current projection (It is not magic you need to check the meta data 
## or ask your instructor) 
## ETRS89 / UTM zone 32N
## the definition of proj4 strings is kind ob obsolet have a look at the links under section zero
epsg_number = 25832

# switch if lasclip is called
lasclip = FALSE
set.seed(1000)

#  area of interest (central MOF)
xmin = 476174.
ymin = 5631386.
xmax = 478217.
ymax = 5632894.

# test area so called "sap flow halfmoon"
xmin = 477500
ymin = 5631730
xmax = 478350
ymax = 5632500

# define variables for the lidR catalog
chunksize = 250
overlap = 25



# get viridris color palette
pal<-mapview::mapviewPalette("mapviewTopoColors")

# 3 - start code 
#-----------------

#---- this part is for clipping only
if (lasclip){
  # Get all *.las files of a the folder you have specified to contain the original las files
  las_files = list.files(envrmt$path_lidar_org, pattern = glob2rx("*.las"), full.names = TRUE)
  
  
  #---- NOTE OTIONALLY you can cut the original big data set  to the smaller extent of the MOF
  # https://www.rdocumentation.org/packages/lidR/versions/1.6.1/topics/lasclip
  core_aoimof<- lidR::lasclipRectangle(lidR::readLAS(las_files[1]), xleft = xmin, ybottom = ymin, xright = xmax, ytop = ymax)
  
  # write the new dataset to the level0 folder and create a corresponding index file (lax)
  lidR::writeLAS(core_aoimof,file.path(envrmt$path_level0,"las_mof.las"))
  rlas::writelax(file.path(envrmt$path_level0, "las_mof.las"))
} #---- OPTIONAL clipping section finished


#---- We assume that you have the "las_mof.las" file in the folder 
#     that is stored in envrmt$path_level0 
# setting up a lidR catalog structure
mof100_ctg <- lidR::readLAScatalog(envrmt$path_level0)
projection(mof100_ctg) <- epsg_number
lidR::opt_chunk_size(mof100_ctg) = chunksize
future::plan(multisession)
lidR::opt_chunk_buffer(mof100_ctg) <- overlap
lidR::opt_output_files(mof100_ctg) <- paste0(envrmt$path_normalized,"/{ID}_norm") # add output filname template
mof100_ctg@output_options$drivers$Raster$param$overwrite <- TRUE

#---- derive DTM DEM and CHM information from an ALS point cloud

# the fastest and simplest algorithm to interpolate a surface is given with p2r()
# the available options are  p2r, dsmtin, pitfree
dsm_p2r_1m = grid_canopy(mof100_ctg, res = 1, algorithm = p2r())
plot(dsm_p2r_1m)              

# now we calculate a digital terrain model by interpolating the ground points 
# and creates a rasterized digital terrain model. The algorithm uses the points 
# classified as "ground" and "water (Classification = 2 and 9 according to LAS file format 
# available algorithms are  knnidw, tin, and kriging
dtm_knnidw_1m <- grid_terrain(mof100_ctg, res=1, algorithm = knnidw(k = 6L, p = 2))
plot(dtm_knnidw_1m)

# we remove the elevation of the surface from the catalog data and create a new catalog
mof100_ctg_chm <- lidR::normalize_height(mof100_ctg, dtm_knnidw_1m)

# if you want to save this catalog  an reread it  you need 
# to uncomment the following lines
#- saveRDS(mof100_ctg_chm,file= file.path(envrmt$path_level1,"mof100_ctg_chm.rds"))
#- mof100_ctg_chm<-readRDS(file.path(envrmt$path_level1,"mof100_ctg_chm.rds"))


# Now create a CHM based on the normalized data and a CHM with the dsmtin() algorithm
# Note the DSM (digital surface model) is now a CHM because we 
# already have normalized the data

# first set a NEW  output name for the catalog
lidR::opt_output_files(mof100_ctg_chm) <- paste0(envrmt$path_normalized,"/{ID}_chm_dsmtin") # add output filname template

# calculate a chm raster with dsmtin()/p2r
chm_dsmtin_1m = grid_canopy(mof100_ctg_chm, res=1.0, dsmtin())

# write it to a tif file
raster::writeRaster(chm_dsmtin_1m,file.path(envrmt$path_data_mof,"chm_dsmtin_1m.tif"),overwrite=TRUE) 

# 4 - visualize 
# -------------------

## standard plot command

plot(raster::raster(file.path(envrmt$path_data_mof,"chm_dsmtin_1m.tif")))

## call mapview with some additional arguments
mapview(raster(file.path(envrmt$path_data_mof,"chm_dsmtin_1m.tif")),
        map.types = "Esri.WorldImagery",  
        legend=TRUE, 
        layer.name = "canopy height model",
        col = pal(256),
        alpha.regions = 0.65)
