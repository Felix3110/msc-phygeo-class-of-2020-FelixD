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

# 0 - load additional packages
#-----------------------------
# for an unique combination of all files in the file list
# google expression: rcran unique combination of vector 
# alternative google expression: expand.grid unique combinations
# have a look at: https://stackoverflow.com/questions/17171148/non-redundant-version-of-expand-grid
add_pkgs = c("gtools")


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
                 "msc-phygeo-class-of-2020-creu/src/000_setup.R"))


# 2 - define variables
#---------------------

## define current projection (It is not magic you need to check the meta data 
## or ask your instructor) 
## ETRS89 / UTM zone 32N
## the definition of proj4 strings is kind ob obsolet have a look at the links under section zero
epsg_number = 25832

set.seed(1000)

# creating a list for files to be deleted if necessary
r_flist = ""

# test area so called "sap flow halfmoon"
xmin = 477500
ymin = 5631730
xmax = 478350

ymax = 5632500


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


for (i in 1:nrow(df)) {
  if (raster(df[i,1])@extent==raster(df[i,2])@extent){ # compare the extent
    cat("fix ",df[i,1]," ", df[i,2],"\n")        # output for information 
    new_raster = stack(df[i,1]) + stack(df[i,2]) - 255  # formula to fix
    cat("write ",paste0(envrmt$path_aerial_org,basename(max(df[i,]))),"\n") # output for information
    writeRaster(new_raster,  paste0(envrmt$path_aerial_org,basename(max(df[i,]))),overwrite=T) # save it
    cat("rename ",paste0(envrmt$path_aerial_org,basename(min(df[i,]))),"\n") # output for information
    r_flist = append(r_flist,paste0(envrmt$path_aerial_org,basename(min(df[i,]))))
    fixed = fixed + 1
  } 
}
cat(no_comb ," combinations checked\n ",fixed," images are fixed\n")
file.remove(r_flist)

#---- now we need to merge these images to on big image
#     again lets have a look what google tell us
#     google expression: merging multiple rasters R cran     
#     https://stackoverflow.com/questions/15876591/merging-multiple-rasters-in-r
#     with regards to Flo Detsch who asked it the right way

# ok lets follow the rabbit

# create a list for the files to be merged
mofimg=list()
# get new filelist
tif_files = list.files(envrmt$path_aerial_org, pattern = glob2rx("*.tif"), full.names = TRUE)
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

raster::writeRaster(sapflow,file.path(envrmt$path_data_mof,"sapflow.tif"),overwrite=TRUE) 
