#installiert envimaR Paket von github
install.packages("devtools")
install.packages("usethis")
install.packages("link2GI")
install.packages("sp")
install.packages("lidR")
install.packages("mapview")
install.packages("rlas")
devtools::install_github("envima/envimaR", force = T)
###AB HIER bis source file path msc laden
library("raster")
library("envimaR")
library("devtools")
library("usethis")
library("link2GI")
library("lidR")
library("mapview")
library("rlas")

rasterOptions(tmpdir = envrmt$path_tmp)
rasterOptions()

otb_path = findOTB()
findOTB()

########-----------------------------------


root_folder = envimaR::alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F:\\BEN\\edu")

project_folders = c("data/",                                 # data folders
                    "data/aerial/org/", "data/lidar/org/", "data/grass/", 
                    "data/data_mof", "data/tmp/", 
                    "run/", "log/",                          # bins and logging
                    "msc-phygeo-class-of-2020-FelixD/src/",   # source code
                    "msc-phygeo-class-of-2020-FelixD/doc/")   # markdown etc. 

# Set libraries
libs = c("link2GI", "raster", "rgdal", "sp", "sf")

# Automatically set root direcory, folder structure and load libraries
envrmt = createEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                    folders = project_folders, 
                    path_prefix = "path_", libs = libs,
                    alt_env_id = "COMPUTERNAME", alt_env_value = "PCRZP",
                    alt_env_root_folder = "F:\\BEN\\edu")

file.path(envrmt$path_src,"000_setup.R")

source(file.path(envrmt$path_src,"000_setup.R"))

source(file.path(root_folder, "msc-phygeo-class-of-2020-FelixD/src/000_setup.R"))


source("C:/Users/Felix/Documents/edu/mpg-envinsys-plygrnd/msc-phygeo-class-of-2020-FelixD/src/000_setup.R")


"~/edu/mpg-envinsys-plygrnd/"
"msc-phygeo-class-of-2020-FelixD"
"/src/000_setup.R"

########

source("https://raw.githubusercontent.com/Jean-Romain/PointCloudViewer/master/sdl.R")  ####pointcloudviewer installieren
devtools::install_github("Jean-Romain/PointCloudViewer")
library("PointCloudViewer")

###Alle Pakete bis hierhin eingeladen (12.11.2020) (Funzt)##########################################################









####folgend die arbeitsumgebung reingeladen


utils::download.file(url="https://github.com/gisma/gismaData/raw/master/uavRst/data/lidR_data.zip",
                     destfile=paste0(envrmt$path_tmp,"chm.zip"))
unzip(paste0(envrmt$path_tmp,"chm.zip"),
      exdir = envrmt$path_tmp,  
      overwrite = TRUE)

las_files = list.files(envrmt$path_tmp,
                       pattern = glob2rx("*.las"),
                       full.names = TRUE)

lidar_file = readLAS(las_files[1])
plot(lidar_file, bg = "green", color = "Z",colorPalette = mvTop(256),backend="pcv")

lidR::plot(las,color = "Z", colorPalette = pal(100),backend = "pcv") #### funzt alles und ist Grundlage für alles weitere









####Creating a Canopy Height Model (CHM) Part 2#######

library("future")


# 1 - source files
#-----------------

#---- source setup file
source(file.path(envimaR::alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd",
                                          alt_env_id = "COMPUTERNAME",
                                          alt_env_value = "PCRZP",
                                          alt_env_root_folder = "F:/BEN/edu"),
                 "src/mpg_course_basic_setup.R"))


# 2 - define variables
#---------------------

## define current projection (It is not magic you need to check the meta data or ask your instructor) 
## ETRS89 / UTM zone 32N
proj4 = "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"


## get viridris color palette
pal<-mapview::mapviewPalette("mapviewTopoColors")

# 3 - start code 
#-----------------


#---- NOTE file size is about 12MB
utils::download.file(url="https://github.com/gisma/gismaData/raw/master/uavRst/data/lidR_data.zip",
                     destfile=paste0(envrmt$path_tmp,"chm.zip"))
unzip(paste0(envrmt$path_tmp,"chm.zip"),
      exdir = envrmt$path_tmp,  
      overwrite = TRUE)

#---- Get all *.las files of a folder into a list
las_files = list.files(envrmt$path_tmp,
                       pattern = glob2rx("*.las"),
                       full.names = TRUE)

#---- create CHM as provided by
# https://github.com/Jean-Romain/lidR/wiki/Rasterizing-perfect-canopy-height-models
las = readLAS(las_files[1])
las = lidR::lasnormalize(las, knnidw())

# reassign the projection
sp::proj4string(las) <- sp::CRS(proj4)

# calculate the chm with the pitfree algorithm
chm = lidR::grid_canopy(las, 0.25, pitfree(c(0,2,5,10,15), c(0,1), subcircle = 0.2))

# write it to tif
raster::writeRaster(chm,file.path(envrmt$path_data_mof,"mof_chm_one_tile.tif"),overwrite=TRUE) 


# 4 - visualize 
-------------------
  
  # call mapview with some additional arguments
  mapview(raster::raster(file.path(envrmt$path_data_mof,"mof_chm_one_tile.tif")),
          legend=TRUE, 
          layer.name = "canopy height model",
          col = pal(256),
          alpha.regions = 0.7)


################Aufgabe 3 noch nicht gelöst


# 0 - load packages
#-----------------------------
library("future")

# 1 - source files
#-----------------
source(file.path(envimaR::alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd",
                                          alt_env_id = "COMPUTERNAME",
                                          alt_env_value = "PCRZP",
                                          alt_env_root_folder = "F:/BEN/edu/mpg-envinsys-plygrnd"),
                 "src/mpg_course_basic_setup.R"))


# 2 - define variables
#---------------------


#  area of interest (central MOF)
xmin<-476174.
ymin<-5631386.
xmax<-478217. 
ymax<-5632894.

# test area so called "sap flow halfmoon"
xmin=477500.0
ymin=5631730.0
xmax=478350.0
ymax=5632500.0

## define variables for the lidR catalog
chunksize = 100
overlap= 10

## define current projection (It is not magic you need to check the meta data or ask your instructor) 
## ETRS89 / UTM zone 32N
proj4 = "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"


# get viridris color palette
pal<-mapview::mapviewPalette("mapviewTopoColors")

# 3 - start code 
#-----------------

# Get all *.las files of a folder into a list
las_files = list.files(envrmt$path_lidar_org, pattern = glob2rx("*.las"), full.names = TRUE)

#---- if you run into memory shortages set up the a lidR catalog 
# we need it for better handling poor available memory 
# check on https://rdrr.io/cran/lidR/man/catalog.html

#---- If not cut the original data to new extent 
# note if using the catalog exchange lidR::readLAS(las_files[1]) with core_aoimof_ctg
# check on https://rdrr.io/cran/lidR/man/lasclip.html
core_aoimof<- lidR::lasclipRectangle(lidR::readLAS(las_files[1]), xleft = xmin, ybottom = ymin, xright = xmax, ytop = ymax)

#---- write the new dataset to the level0 folder and create a corresponding index file (lax)
lidR::writeLAS(core_aoimof,file.path(envrmt$path_level0,"las_mof.las"))
rlas::writelax(file.path(envrmt$path_level0, "las_mof.las"))

#---- setting up lidR catalog
mof100_ctg <- lidR::readLAScatalog(envrmt$path_level0)
sp::proj4string(mof100_ctg) <- sp::CRS(proj4)
future::plan(multisession, workers = 2L)
lidR::set_lidr_threads(4)
lidR::opt_chunk_size(mof100_ctg) = chunksize
lidR::opt_chunk_buffer(mof100_ctg) <- overlap
lidR::opt_output_files(mof100_ctg) <- paste0(envrmt$path_normalized,"/{ID}_norm") # add output filname template
mof100_ctg@output_options$drivers$Raster$param$overwrite <- TRUE

#---- calculate chm following example 1 from the help
# Remove the topography from a point cloud 
dtm <- lidR::grid_terrain(mof100_ctg, 1, lidR::kriging(k = 10L))
mof100_ctg_norm <- lidR::lasnormalize(mof100_ctg, dtm)

# if you want to save this catalog  an reread it  you need to uncomment the following lines
# saveRDS(mof100_ctg_norm,file= file.path(envrmt$path_level1,"mof100_ctg_norm.rds"))
# mof100_ctg_norm<-readRDS(file.path(envrmt$path_level1,"mof100_ctg_norm.rds"))


# Create a CHM based on the normalized data and a DSM by pitfree()
# if the error  "filename exists; use overwrite=TRUE" occures navigate to the 
# paste0(envrmt$path_normalized,"/{ID}_norm") and delete the tif files

lidR::opt_output_files(mof100_ctg_norm) <- paste0(envrmt$path_normalized,"/{ID}_chm") # add output filname template
chm_all = grid_canopy(mof100_ctg_norm, 1.0, dsmtin())

# alternative using pitfree
#chm_all = grid_canopy(mof100_ctg_norm_kri, 1.0, pitfree(c(0,2,5,10,15), c(0,1)))

# write it to tiff
raster::writeRaster(chm_all,file.path(envrmt$path_data_mof,"mof_chm_all.tif"),overwrite=TRUE) 

# 4 - visualize 
-------------------
  
  ## set mapview raster options to full resolution
  
  ## visualize the catalogs
  mapview(mof100_ctg) + mapview(mof_norm, zcol= "Max.Z" )

## call mapview with some additional arguments
mapview(raster::raster(file.path(envrmt$path_data_mof,"mof_chm_all.tif")),
        legend=TRUE, 
        layer.name = "canopy height model",
        col = pal(256),
        alpha.regions = 0.65)




########

#------------------------------------------------------------------------------
# Type: control script 
# Name: 10_preprocess_RGB.R
# Author: Chris Reudenbach, creuden@gmail.com
# Description:  - fixes the white areas of airborne RGB images
#               - merge all image to one image
#               - clip the image to the sapflow half moon
#              
#              
# Data: regular authority provided airborne RGB imagery 
# Output: merged, clipped and corrected image of AOI
# Copyright: Chris Reudenbach, Thomas Nauss 2017,2020, GPL (>= 3)
# git clone https://github.com/GeoMOER-Students-Space/msc-phygeo-class-of-2020-creu.git
#------------------------------------------------------------------------------

## clean your environment
rm(list=ls()) 

# 0 - load packages
#-----------------------------
# for an unique combination of all files in the file list
# google expression: rcran unique combination of vector 
# alternative google expression: expand.grid unique combinations
# have a look at: https://stackoverflow.com/questions/17171148/non-redundant-version-of-expand-grid
install.packages("gtools")
library("gtools")


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


# for an unique combination of all files in the file list
# google expression: rcran unique combination of vector 
# alternative google expression: expand.grid unique combinations
# have a look at: https://rdrr.io/cran/gimme/src/R/expand.grid.unique.R
expand.grid.unique <- function(x, y, incl.eq = FALSE){
  g <- function(i){
    z <- setdiff(y, x[seq_len(i - incl.eq)])
    if(length(z)) cbind(x[i], z, deparse.level = 0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}



# 3 - start code 
#-----------------

#---- this part is for clipping only

# Get all *.las files of a the folder you have specified to contain the original las files
tif_files = list.files(envrmt$path_aerial_org, pattern = glob2rx("*.tif"), full.names = TRUE)

# we create a unique combination of all files for a brute force comparision of the extents
# first possibility using the sourced function from the gimme library
df = expand.grid.unique(tif_files,tif_files)
# second possibility using the combinations() function from the gtools package
#----- NOTE it is not necessary to write the code by yourself 
#      you need to know (1) What do you want to do (2) What do you need (3) express it to google
df = combinations(n = length(tif_files), r = 2, v = tif_files, repeats.allowed = FALSE)

#---- The idea is to check for each unique pair of files if they have the same extent
#     if so we fix the white area using the formula image A + image B - 255
#     additionally we write this file to the name without the "_1" extension and 
#     rename the "_1" file in "_1~" what t least means on unix flavoured systems 
#     that it is handled as a backup file

no_comb = nrow(df)
fixed = 0
# for loop for each element of the data frame (nrow())
for (i in 2:nrow(df)) {
  if (raster(df[i,1])@extent==raster(df[i,2])@extent){ # compare the extent
    cat("fix ",df[i,1]," ", df[i,2],"\n")        # output for information 
    new_raster = stack(df[i,1]) + stack(df[i,2]) - 255  # formula to fix
    cat("write ",paste0(envrmt$path_aerial_org, "/", basename(max(df[i,]))),"\n") # output for information
    writeRaster(new_raster,  paste0(envrmt$path_aerial_org,"/", basename(max(df[i,]))),overwrite=T) # save it
    cat("rename ",paste0(envrmt$path_aerial_org,"/", basename(min(df[i,]))),"\n") # output for information
    file.rename(paste0(envrmt$path_aerial_org,"/", basename(min(df[i,]))),paste0(envrmt$path_aerial_org,"/", basename(min(df[i,])),"~")) # rename
    fixed = fixed + 1
  } 
  # due to renaming we have to rebuild the list        
  tif_files = list.files(envrmt$path_aerial_org, pattern = glob2rx("*.tif"), full.names = TRUE)
  df = expand.grid.unique(tif_files,tif_files)                
}
cat(no_comb ," combinations checked\n ",fixed," images are fixed")


#---- now we need to merge these images to on big image
#     again lets have a look what google tell us
#     google expression: merging multiple rasters R cran     
#     https://stackoverflow.com/questions/15876591/merging-multiple-rasters-in-r
#     with regards to Flo Detsch who asked it the right way

# ok lets follow the rabbit
# creating a list
mofimg=list()
# stacking all files and put them in the list object
cat("stack files...\n")
for (f in 1:length(tif_files)){
  mofimg[[f]] = stack(tif_files[f])
}

# setting the parameters
mofimg$tolerance = 1
mofimg$filename  = paste0(envrmt$path_aerial_org,"merged_mof.tif")
mofimg$overwrite = TRUE
cat("merge files - this will take a while \n")
merged_mof = do.call(raster::merge, mofimg)

# cropping it using the ?crop or ?clip and decide  
# or google for something like "crop clip raster images R cran"
# as a result you need a vector file or an extent

# defining the extent object using the above defined params 
# for more info ?extent
cat("crop AOI\n")
ext <- extent(xmin,xmax,ymin,ymax)

# cropping it
sapflow  =  crop(merged_mof, ext)	

# 4 - visualize 
# -------------------

plotRGB(merged_mof)

plotRGB(sapflow)


