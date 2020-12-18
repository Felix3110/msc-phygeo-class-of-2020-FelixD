#------------------------------------------------------------------------------
# Type: control script 
# Name: 20_segmentation_Catalog.R
# Author: Chris Reudenbach, creuden@gmail.com
# Description:  calculates a comprehensive set of tree segmentations based on the CHM data set
# Data: CHM raster file as derived by 10_CHM_Catalog.R 
# Output: Segmentation layers
# Copyright: Chris Reudenbach, Thomas Nauss 2017,2020, GPL (>= 3)
#------------------------------------------------------------------------------

## clean your environment
rm(list=ls()) 

# 0 - load packages
#-----------------------------


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
  las_files = list.files(envrmt$path_level0,"las_mof.las", pattern = glob2rx("*.las"), full.names = TRUE)
  
  
  
  #---- We assume that you have the "las_mof.las" file in the folder 
  #     that is stored in envrmt$path_level0 
  # setting up a lidR catalog structure
  mof100_ctg <- lidR::readLAScatalog(envrmt$path_level0)
  projection(mof100_ctg) <- epsg_number
  lidR::opt_chunk_size(mof100_ctg) = chunksize
  future::plan(multisession)
  lidR::opt_chunk_buffer(mof100_ctg) <- overlap
  lidR::opt_output_files(mof100_ctg) <- paste0(envrmt$path_level0,"{ID}_norm") # add output filname template
  mof100_ctg@output_options$drivers$Raster$param$overwrite <- TRUE
  
  
  # 4 - visualize 
  # -------------------
  
  ## standard plot command
  plot(raster::raster(file.path(envrmt$path_data_mof,"chm_dsmtin_1m.tif")))
  mapview(mof100_ctg) + mapview(mof100_ctg, zcol= "Max.Z" )
  