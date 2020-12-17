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
