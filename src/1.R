#installiert envimaR Paket von github
install.packages("devtools")
install.packages("usethis")
install.packages("link2GI")
devtools::install_github("envima/envimaR", force = T)
library("raster")
library("envimaR")
library("devtools")
library("usethis")
library("link2GI")
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
                    "msc-phygeo-class-of-2020-FelixD./src/",   # source code
                    "msc-phygeo-class-of-2020-FelixD./doc/")   # markdown etc. 

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

source(file.path(root_folder, "msc-phygeo-class-of-2020-FelixD./000_setup.R"))


source("C:/Users/Felix/Documents/edu/mpg-envinsys-plygrnd/msc-phygeo-class-of-2020-FelixD/src/000_setup.R")


"~/edu/mpg-envinsys-plygrnd/"
"msc-phygeo-class-of-2020-FelixD"
"/src/000_setup.R"

###Alle Pakete bis hierhin eingeladen (12.11.2020) (Funzt)

