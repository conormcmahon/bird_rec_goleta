

# ********************************************************
# ***************** Load Data and Libraries ******************

library(dplyr)
library(raster)
library(here)

# Load raster data, corresponding to
#   Arroyo Quemado (Baron Ranch)
files_aq <- list.files(path="D:/birdrec/remote_sensing_data/arroyo_quemado/LiDAR/pcd/output/vegetation_histograms/", pattern="*histogram.tif")
# Devereux Creek (and others, on Coal Oil Point Reserve, North Campus Open Space, and Ellwood Mesa)
files_copr <- list.files(path="D:/birdrec/remote_sensing_data/devereux_ellwoood/LiDAR/pcd/output/vegetation_histograms/", pattern="*histogram.tif")
#   Figueroa Creek (Sedgwick Ranch)
files_sr <- list.files(path="D:/SERDP/SHIFT/LiDAR_Sedgwick/pcd/output/vegetation_histograms/", pattern="*histogram.tif")
#   Various Creeks (Santa Ynez Mountains - Front Range, or South Slope)
files_sym <- list.files(path="D:/SERDP/GEE_Classifier/LiDAR/pcd/output/vegetation_histograms/", pattern="*histogram.tif")
#   Atascadero Creek 
files_ac <- list.files(path="D:/birdrec/remote_sensing_data/atascadero/LiDAR/pcd/output/vegetation_histograms/", pattern="*histogram.tif")
#   Maria Ygnacio Creek 
files_myc <- list.files(path="D:/birdrec/remote_sensing_data/maria_ygnacio/LiDAR/pcd/output/vegetation_histograms/", pattern="*histogram.tif")
#   Cold Spring Creek
files_csc <- list.files(path="D:/SERDP/Mission_Creek/LiDAR_cold_spring/pcd/output/vegetation_histograms/", pattern="*histogram.tif")




# ********************************************************
# ***************** Mosaicking Function ******************

# Add a new raster stack (from a file path) to an existing raster mosaic
#   At overlapping points, retains the maximum 
addImageToMosaic <- function(new_image_path, existing_image)
{
  print(paste("  Loading a new image from ",
              new_image_path), sep="")
  # If existing image is empty, just return the new image
  if(!hasValues(existing_image))
    return(stack(new_image_path))
  # Else, mosaic the two together
  raster::mosaic(stack(new_image_path), 
                 existing_image, 
                 fun=max)
}


# ********************************************************
# ***************** Generate Arroyo Quemado Mosaic *****************

print(paste("Starting to construct raster for Arroyo Quemado data, from ",
            length(files_aq), 
            " input rasters.", sep=""))
mosaic_aq <- raster()
# There's probably a better way to do this without a for loop but I am too tired to think of it 
for(file in files_aq)
{
  mosaic_aq <- addImageToMosaic(paste("D:/birdrec/remote_sensing_data/arroyo_quemado/LiDAR/pcd/output/vegetation_histograms", file, sep="/"),
                                  mosaic_aq)
}
writeRaster(mosaic_aq, "D:/birdrec/remote_sensing_data/arroyo_quemado/LiDAR/pcd/output/vegetation_histograms/vegetation_histogram_mosaic.tif", overwrite=TRUE)




# ********************************************************
# ***************** Generate COPR / Ellwood Mesa / NCOS Mosaic *****************

print(paste("Starting to construct raster for COPR / Ellwood Mesa / NCOS data, from ",
            length(files_copr), 
            " input rasters.", sep=""))
mosaic_copr <- raster()
# There's probably a better way to do this without a for loop but I am too tired to think of it 
for(file in files_copr)
{
  mosaic_copr <- addImageToMosaic(paste("D:/birdrec/remote_sensing_data/devereux_ellwoood/LiDAR/pcd/output/vegetation_histograms", file, sep="/"),
                                  mosaic_copr)
}
writeRaster(mosaic_copr, "D:/birdrec/remote_sensing_data/devereux_ellwoood/LiDAR/pcd/output/vegetation_histograms/vegetation_histogram_mosaic.tif", overwrite=TRUE)






# ********************************************************
# ***************** Generate Sedgwick Ranch Mosaic *****************

print(paste("Starting to construct raster for sr data, from ",
            length(files_sr), 
            " input rasters.", sep=""))
mosaic_sr <- raster()
# There's probably a better way to do this without a for loop but I am too tired to think of it 
for(file in files_sr)
{
  mosaic_sr <- addImageToMosaic(paste("D:/SERDP/SHIFT/LiDAR_Sedgwick/pcd/output/vegetation_histograms", file, sep="/"),
                                  mosaic_sr)
}
writeRaster(mosaic_sr, "D:/SERDP/SHIFT/LiDAR_Sedgwick/pcd/output/vegetation_histograms/vegetation_histogram_mosaic.tif", overwrite=TRUE)




# ********************************************************
# ***************** Generate Santa Ynez Mountain Front Range Mosaic *****************

print(paste("Starting to construct raster for sym data, from ",
            length(files_sym), 
            " input rasters.", sep=""))
mosaic_sym <- raster()
# There's probably a better way to do this without a for loop but I am too tired to think of it 
for(file in files_sym)
{
  mosaic_sym <- addImageToMosaic(paste("D:/SERDP/GEE_Classifier/LiDAR/pcd/output/vegetation_histograms", file, sep="/"),
                                  mosaic_sym)
}
writeRaster(mosaic_sym, "D:/SERDP/GEE_Classifier/LiDAR/pcd/output/vegetation_histograms/vegetation_histogram_mosaic.tif", overwrite=TRUE)









# ********************************************************
# ***************** Generate Atascadero Creek Mosaic *****************

print(paste("Starting to construct raster for ac data, from ",
            length(files_ac), 
            " input rasters.", sep=""))
mosaic_ac <- raster()
# There's probably a better way to do this without a for loop but I am too tired to think of it 
for(file in files_ac)
{
  mosaic_ac <- addImageToMosaic(paste("D:/birdrec/remote_sensing_data/atascadero/LiDAR/pcd/output/vegetation_histograms", file, sep="/"),
                                mosaic_ac)
}
writeRaster(mosaic_ac, "D:/birdrec/remote_sensing_data/atascadero/LiDAR/pcd/output/vegetation_histograms/vegetation_histogram_mosaic.tif", overwrite=TRUE)




# ********************************************************
# ***************** Generate Maria Ygnacio Creek Mosaic *****************

print(paste("Starting to construct raster for myc data, from ",
            length(files_myc), 
            " input rasters.", sep=""))
mosaic_myc <- raster()
# There's probably a better way to do this without a for loop but I am too tired to think of it 
for(file in files_myc)
{
  mosaic_myc <- addImageToMosaic(paste("D:/birdrec/remote_sensing_data/maria_ygnacio/LiDAR/pcd/output/vegetation_histograms", file, sep="/"),
                                 mosaic_myc)
}
writeRaster(mosaic_myc, "D:/birdrec/remote_sensing_data/maria_ygnacio/LiDAR/pcd/output/vegetation_histograms/vegetation_histogram_mosaic.tif", overwrite=TRUE)




# ********************************************************
# ***************** Generate Cold Spring Creek Mosaic *****************

print(paste("Starting to construct raster for csc data, from ",
            length(files_csc), 
            " input rasters.", sep=""))
mosaic_csc <- raster()
# There's probably a better way to do this without a for loop but I am too tired to think of it 
for(file in files_csc)
{
  mosaic_csc <- addImageToMosaic(paste("D:/SERDP/Mission_Creek/LiDAR_cold_spring/pcd/output/vegetation_histograms", file, sep="/"),
                                 mosaic_csc)
}
writeRaster(mosaic_csc, "D:/SERDP/Mission_Creek/LiDAR_cold_spring/pcd/output/vegetation_histograms/vegetation_histogram_mosaic.tif", overwrite=TRUE)
