1
#######
devtools::install_github("envima/envimaR", force = T)
devtools::install_github("gisma/uavRst")
devtools::install_github("r-spatial/link2GI", force = T)
library(envimaR)
######
packagesToLoad = c("lidR", "link2GI", "mapview", "raster", "rgdal", "rlas", "sp",  "sf")
#####
mvTop<-mapview::mapviewPalette("mapviewTopoColors")
mvSpec<-mapviewTopoColors<-mapview::mapviewPalette("mapviewSpectralColors")
#####
pal<-mapview::mapviewPalette("mapviewTopoColors")
###################################################################################################


rootDir = envimaR::alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd/",
                                   alt_env_id = "COMPUTERNAME",
                                   alt_env_value = "PCRZP",
                                   alt_env_root_folder = "F:/BEN/edu")

#####
projectDirList   = c("data/",                # datafolders for all kind of date
                     "data/auxdata/",        # the following used by scripts however
                     "data/aerial/",     # you may add whatever you like                     
                     "data/aerial/org/",     # you may add whatever you like
                     "data/lidar/org/",
                     "data/lidar/",
                     "data/grass/",
                     "data/lidar/level0/",                     
                     "data/lidar/level1/",
                     "data/lidar/level1/normalized",
                     "data/lidar/level1/ID",
                     "data/lidar/level2/",
                     "data/lidar/level0/all/",
                     "data/data_mof", 
                     "data/tmp/",
                     "run/",                # temporary data storage
                     "log/",                # logging
                     "msc-phygeo-class-of-2020-FelixD/src/",                # scripts
                     "msc-phygeo-class-of-2020-FelixD/doc/")                # documentation markdown etc.

#####
# setup of root directory, folder structure and loading libraries
# returns "envrmt" list which contains the folder structure as short cuts
envrmt = envimaR::createEnvi(root_folder = rootDir,
                             folders = projectDirList,
                             path_prefix = "path_",
                             libs = packagesToLoad,
                             alt_env_id = "COMPUTERNAME",
                             alt_env_value = "PCRZP",
                             alt_env_root_folder = "F:/BEN/edu")

### set raster temp path
raster::rasterOptions(tmpdir = envrmt$path_tmp)
return(envrmt)

###########################

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
##################################################
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

######OTB
otb_path = findOTB()
findOTB()
