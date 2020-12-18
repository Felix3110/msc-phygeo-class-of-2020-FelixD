#------------------------------------------------------------------------------
# Type: control script 
# Name: 20_calculate_synthetic_bands.R
# Author: Chris Reudenbach, creuden@gmail.com
# Description:  - calculate on base of the original images
#               - RGB Indices
#               - structural images
#               - statistical derivations
#              
#              
# Data: corrected RGB image of AOI 
# Output: comprehensive image stack of useful synthetic bands 
# Copyright: Chris Reudenbach, Thomas Nauss 2017,2020, GPL (>= 3)
# git clone https://github.com/GeoMOER-Students-Space/msc-phygeo-class-of-2020-creu.git
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
source(file.path(envimaR::alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd/",
                                          alt_env_id = "COMPUTERNAME",
                                          alt_env_value = "PCRZP",
                                          alt_env_root_folder = "F:/BEN/edu/mpg-envinsys-plygrnd"),
                 "src/000_setup.R"))


# 2 - define variables
#---------------------

## define current projection (It is not magic you need to check the meta data 
## or ask your instructor) 
## ETRS89 / UTM zone 32N
## the definition of proj4 strings is kind ob obsolet have a look at the links under section zero
epsg_number = 25832

set.seed(1000)



# test area so called "sap flow halfmoon"
xmin = 477500
ymin = 5631730
xmax = 478350
ymax = 5632500


# get viridris color palette
pal<-mapview::mapviewPalette("mapviewTopoColors")






# 3 - start code 
#-----------------

#---- this part is for clipping only

# Get all *.las files of a the folder you have specified to contain the original las files
tif_files = list.files(envrmt$path_aerial_org, pattern = glob2rx("*.tif"), full.names = TRUE)




# 4 - visualize 
# -------------------