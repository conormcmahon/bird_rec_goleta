
# Code to sample remote sensing products for vegetation: 
#   Phenology (Sentinel-2 processed in GEE to 12-month greenness timeseries)
#   Structure (LiDAR height histograms)

library(raster)
library(tidyverse)
library(janitor)
library(here)
library(sf)
library(paletteer)
library(rasterKernelEstimates)
library(caret)

# Load phenology rasters
datastack_aq <- brick("D:/birdrec/remote_sensing_data/gee_results/baron_datastack_2021.tif")
datastack_copr <- brick("D:/birdrec/remote_sensing_data/gee_results/copr_ellwood_datastack_2021.tif")
datastack_sr <- brick("D:/birdrec/remote_sensing_data/gee_results/sedgwick_datastack_2021.tif")
datastack_sym <- brick("D:/SERDP/GEE_Classifier/Dar_manual_delineations/sym_fr_datastack_2021.tif")

# Load height histogram rasters
#   NOTE - annoyingly, I did not output the LiDAR rasters in the same CRS as the GEE rasters... so we'll reproject the LiDAR data here
lidar_aq <- projectRaster(brick("D:/birdrec/remote_sensing_data/arroyo_quemado/LiDAR/pcd/output/vegetation_histograms/vegetation_histogram_mosaic.tif"), datastack_aq)
lidar_copr <- projectRaster(brick("D:/birdrec/remote_sensing_data/devereux_ellwoood/LiDAR/pcd/output/vegetation_histograms/vegetation_histogram_mosaic.tif"), datastack_copr)
lidar_sr <- projectRaster(brick("D:/SERDP/SHIFT/LiDAR_Sedgwick/pcd/output/vegetation_histograms/vegetation_histogram_mosaic.tif"), datastack_sr)
lidar_sym <- projectRaster(brick("D:/SERDP/GEE_Classifier/LiDAR/pcd/output/vegetation_histograms/vegetation_histogram_mosaic.tif"), datastack_sym)
lidar_ac <- projectRaster(brick("D:/birdrec/remote_sensing_data/atascadero/LiDAR/pcd/output/vegetation_histograms/vegetation_histogram_mosaic.tif"), datastack_sym)
lidar_myc <- projectRaster(brick("D:/birdrec/remote_sensing_data/maria_ygnacio/LiDAR/pcd/output/vegetation_histograms/vegetation_histogram_mosaic.tif"), datastack_sym)
lidar_ac_myc <- raster::mosaic(lidar_ac, lidar_myc, fun=max)

# Generate simple three-band visualization for LiDAR, binning 5-10m, 10-15m, and 15+m together
lidar_sym_low <- sum(lidar_sym[[8:12]])
lidar_sym_mid <- sum(lidar_sym[[13:17]])
lidar_sym_high <- sum(lidar_sym[[18:32]])
lidar_sym_threeband <- stack(lidar_sym_low,lidar_sym_mid,lidar_sym_high)
lidar_sym_threeband[lidar_sym_threeband] 
# Filter out powerlines:
filterPowerLines <- function(input_raster, kernel_halfwidth, percentile, height_threshold)
{
  maxHeight <- calc(input_raster[[3:nbands(input_raster)]], function(hist_values){
    hist_nonzero <- hist_values > 0
    highest_bin <- max(hist_nonzero * (1:(nbands(input_raster)-2)))
  })
  # Weights matrix is for 24 nearest neighbors
  kernel_width <- (kernel_halfwidth*2+1)
  weights <- matrix(1,nrow=kernel_width, ncol=kernel_width)
  weights[kernel_halfwidth+1,kernel_halfwidth+1] <- 0
  run.time <- proc.time()
  kernel_percentiles <- rasterLocalQuantiles(maxHeight,weights,q=percentile)
  
  powerline_mask <- (maxHeight - kernel_percentiles) > height_threshold
  
  return(mask(input_raster, powerline_mask, maskvalue=1, updatevalue=NA))
}


# Load points shapefile
field_points <- st_read("D:/birdrec/field_data_goleta_2022/goleta_qfield_final/birdrec_site_goleta_2022.gpkg") %>% 
  mutate(creek_name = substring(site_name,1,(nchar(as.character(site_name))-3)))
# Evenly sampled points along Mission and Rattlesnake Creeks
even_points_mission <- st_read("D:/SERDP/Mission_Creek/field_data/mission_chained.gpkg", fid_column_name = "FID") %>%
  mutate(dem = dem_mosaic_sym,
         chm = chm_mosaic_sym,
         slope = slope_mosaic_sym,
         dist_up_creek = as.numeric(as.character(FID))*30) %>%
  select(dem, chm, slope, dist_up_creek, geom)
even_points_mission$creek <- "mission"
even_points_rattlesnake <- st_read("D:/SERDP/Mission_Creek/field_data/rattlesnake_chained.gpkg", fid_column_name = "FID") %>%
  mutate(dem = dem_mosaic_sym,
         chm = chm_mosaic_sym,
         slope = slope_mosaic_sym,
         dist_up_creek = as.numeric(as.character(FID))*30 + 28*30) %>%
  select(dem, chm, slope, dist_up_creek, geom)
even_points_rattlesnake$creek <- "rattlesnake"

# Load locations of pools for Mission and Rattlesnake
pools_rattlesnake_winter <- st_read("D:/SERDP/Mission_Creek/field_data/paul_pools_maps/MCRC_water_winter_2022/RC_pools_water.shp")
pools_rattlesnake_summer <- st_read("D:/SERDP/Mission_Creek/field_data/paul_pools_maps/MCRC_water_summer_2022/RC_pools_summer.shp")
pools_rattlesnake_dry <- st_read("D:/SERDP/Mission_Creek/field_data/paul_pools_maps/MCRC_water_winter_2022/RC_pools_all.shp")
pools_mission_winter <- st_read("D:/SERDP/Mission_Creek/field_data/paul_pools_maps/MCRC_water_winter_2022/MC_pools_water.shp")
pools_mission_summer <- st_read("D:/SERDP/Mission_Creek/field_data/paul_pools_maps/MCRC_water_summer_2022/MC_pools_summer.shp")
pools_mission_dry <- st_read("D:/SERDP/Mission_Creek/field_data/paul_pools_maps/MCRC_water_winter_2022/MC_pools_all.shp")

# For each point in the evenly sampled point set, get the distance to the nearest wet pool in summer, winter, and to the nearest pool overall
even_points_rattlesnake$winter_water_dist <- apply(st_distance(even_points_rattlesnake, st_transform(pools_rattlesnake_winter, st_crs(even_points_rattlesnake))), MARGIN=1, min)
even_points_rattlesnake$summer_water_dist <- apply(st_distance(even_points_rattlesnake, st_transform(pools_rattlesnake_summer, st_crs(even_points_rattlesnake))), MARGIN=1, min)
even_points_rattlesnake$any_pool_dist <- apply(st_distance(even_points_rattlesnake, st_transform(pools_rattlesnake_dry, st_crs(even_points_rattlesnake))), MARGIN=1, min)
even_points_mission$winter_water_dist <- apply(st_distance(even_points_mission, st_transform(pools_mission_winter, st_crs(even_points_mission))), MARGIN=1, min)
even_points_mission$summer_water_dist <- apply(st_distance(even_points_mission, st_transform(pools_mission_summer, st_crs(even_points_mission))), MARGIN=1, min)
even_points_mission$any_pool_dist <- apply(st_distance(even_points_mission, st_transform(pools_mission_dry, st_crs(even_points_mission))), MARGIN=1, min)

even_points <- rbind(even_points_mission, even_points_rattlesnake)
even_points_transformed <- st_cast(st_transform(even_points,crs(lidar_sym)),"POINT")
dist_thresh <- 5
even_points_transformed$summer_dry <- even_points_transformed$any_pool_dist < (even_points_transformed$summer_water_dist-dist_thresh)
even_points_transformed$winter_dry <- even_points_transformed$any_pool_dist < (even_points_transformed$winter_water_dist-dist_thresh)

geologic_boundaries <- data.frame(
  lower_formation = c("Alluvium","Sespe","Coldwater","Cozydell","Matilija"),
  upper_formation = c("Sespe","Coldwater","Cozydell","Matilija","Juncal"),
  elevation_mission = c(180, 275, 550, 700, 930),
  elevation_rattlesnake = c(180, 319, 500, 614, 830),
  position_mission = c(1800, 3570, 6600, 7680, 9000),
  position_rattlesnake = c(1860, 4560, 6330, 7170, 7980),
  ghost_positioning = c(-1,-1,-1,-1,-1)
)

geologic_boundaries_v2 <- data.frame(
  creek = c(rep("Mission",5), rep("Rattlesnake", 5)),
  lower_formation = c("Alluvium","Sespe","Coldwater","Cozydell","Matilija","Alluvium","Sespe","Coldwater","Cozydell","Matilija"),
  upper_formation = c("Sespe","Coldwater","Cozydell","Matilija","Juncal","Sespe","Coldwater","Cozydell","Matilija","Juncal"),
  elevation = c(180, 275, 550, 700, 930, 180, 319, 500, 614, 830),
  position = c(1800, 3570, 6600, 7680, 9000, 1860, 4560, 6330, 7170, 7980),
  ghost_positioning = c(-1,-1,-1,-1,-1)
)

mission_tunnel <- data.frame(
  elevation = 365,
  position = 4860,
  name = "Mission Tunnel",
  ghost_positioning = c(-1)
)
# Conservative boundaries
# summer_water_boundaries <- data.frame(
#   creek = c("Mission","Mission","Rattlesnake","Rattlesnake"),
#   elevation = c(245, 341, 346, 496),
#   position = c(3060, 4560, 4950, 6300)
# )
# Liberal boundaries (including more of reach with patchy surface water)
summer_water_boundaries <- data.frame(
  
  creek = c("Mission","Mission","Rattlesnake","Rattlesnake"),
  elevation = c(230, 341, 294, 496),
  position = c(2820, 4560, 4290, 6300)
)
winter_water_boundaries <- data.frame(
  creek = c("Mission","Mission","Mission","Mission","Rattlesnake","Rattlesnake"),
  elevation = c(202, 339, 511, 830, 269, 466),
  position = c(2250, 4560, 6300, 8040, 870, 6780)
)

# Add information to each point on which formation it overlays 
getFormationAtPoint <- function(index)
{
  target_elev <- even_points_transformed[index,]$dem
  if(even_points_transformed[index,]$creek == "rattlesnake")
    return(as.character(geologic_boundaries[min(which(geologic_boundaries$elevation_rattlesnake > target_elev)),]$lower_formation))
  else
    return(as.character(geologic_boundaries[min(which(geologic_boundaries$elevation_mission > target_elev)),]$lower_formation))
}
even_points_transformed$formation <- unlist(lapply(1:nrow(even_points_transformed), addFormationInfo))

st_write(even_points_transformed, "D:/SERDP/Mission_Creek/field_data/creek_geomorphology_hydrology.gpkg", append=FALSE)

# Check if a value is exactly TRUE (returns 0 if value is NA)
exactlyTrue <- function(x) {
  if(is.na(x)) 
    return(FALSE)
  else
    return(x == 1)
}

# Load densiometer measurements
densiometer_measurements <- st_read("D:/SERDP/Mission_Creek/field_data/2022_09_19_mission/hydro_points_mc_rc.gpkg", fid_column_name = "FID") %>%
  mutate(FID = as.numeric(as.character(FID)))
densiometer_measurements <- densiometer_measurements[((is.na(densiometer_measurements$cover_upst)+is.na(densiometer_measurements$cover_down)+is.na(densiometer_measurements$cover_righ)+is.na(densiometer_measurements$cover_left))) != 4,] %>%
  drop_na(cover_upst, cover_down, cover_righ, cover_left) %>%
  mutate(overall_cover = (floor(cover_upst)+floor(cover_down)+floor(cover_righ)+floor(cover_left))/37/4)  %>%
  mutate(decid_upst = cover_upst - floor(cover_upst),
         decid_down = cover_down - floor(cover_down),
         decid_righ = cover_righ - floor(cover_righ),
         decid_left = cover_left - floor(cover_left)) %>%
  mutate(decid_cover = (decid_down+decid_upst+decid_righ+decid_left)*100/37/4)
densiometer_measurements$summer_water_dist <- apply(st_distance(densiometer_measurements, even_points_transformed %>% filter(!summer_dry)), MARGIN=1, min)
densiometer_measurements$winter_water_dist <- apply(st_distance(densiometer_measurements, even_points_transformed %>% filter(!winter_dry)), MARGIN=1, min)
densiometer_measurements <- st_cast(st_transform(densiometer_measurements,crs(lidar_sym)),"POINT")
densiometer_measurements[unlist(lapply(densiometer_measurements$surf_water, exactlyTrue)),]$winter_water_dist <- 0
densiometer_measurements[unlist(lapply(densiometer_measurements$surf_water, exactlyTrue)),]$summer_water_dist <- 0
densiometer_measurements_decid_est <- densiometer_measurements %>% filter(FID >= 65)

st_write(densiometer_measurements, "D:/SERDP/Mission_Creek/field_data/2022_09_19_mission/cover_estimates.gpkg")

t.test((densiometer_measurements_decid_est %>% filter(summer_water_dist < 5))$decid_cover, 
       (densiometer_measurements_decid_est %>% filter(summer_water_dist >= 5))$decid_cover)
t.test((densiometer_measurements %>% filter(summer_water_dist < 5))$overall_cover, 
       (densiometer_measurements %>% filter(summer_water_dist >= 5))$overall_cover)

# Summary statistics on overall cover
min(densiometer_measurements$overall_cover)
max(densiometer_measurements$overall_cover)
mean(densiometer_measurements$overall_cover)
sd(densiometer_measurements$overall_cover)
# Summary statistics on deciduous cover
min(densiometer_measurements_decid_est$decid_cover)
max(densiometer_measurements_decid_est$decid_cover)
mean(densiometer_measurements_decid_est$decid_cover)
sd(densiometer_measurements_decid_est$decid_cover)
# Overall Cover - plot of variation
densiometer_overall_hist <- ggplot(densiometer_measurements) + 
  geom_histogram(aes(x=overall_cover, y=stat(density)), binwidth=0.05) + 
  xlab("Overall Riparian Canopy Cover in Densiometer Measurements") + 
  ylab("Percent of Sites") + 
  ggtitle("Variation in Canopy Cover Across Densiometer Sites")
ggsave("D:/SERDP/Mission_Creek/output_imagery/canopy_cover_overall_hist.png", densiometer_overall_hist, width=5, height=2.5, device="png")
# Deciduous Cover - plot of variation
densiometer_overall_hist <- ggplot(densiometer_measurements_decid_est) + 
  geom_histogram(aes(x=decid_cover, y=stat(density)), binwidth=0.05) + 
  xlab("Deciduous Riparian Canopy Cover in Densiometer Measurements") + 
  ylab("Percent of Sites") + 
  ggtitle("Variation in Deciduous Cover Across Densiometer Sites")
ggsave("D:/SERDP/Mission_Creek/output_imagery/canopy_cover_decid_hist.png", densiometer_overall_hist, width=5, height=2.5, device="png")


# Load field vegetation maps 
veg_polygons <- st_read("D:/SERDP/Mission_Creek/field_data/2022_09_19_mission/veg_polygons.shp") %>%
  mutate(class = as.character(level_1))
veg_polygons[veg_polygons$class == "Rip",]$class = "Deciduous"
veg_polygons[veg_polygons$class == "Riparian",]$class = "Deciduous"
veg_polygons[veg_polygons$class == "EvergreenTree",]$class = "Evergreen"


# Some basic visualizations of creek properties:
#  DEM vs Creek Position
creek_slope_plot <- ggplot() + 
  geom_point(data=even_points_transformed, aes(x=dist_up_creek, y=dem, group=creek, col=creek), alpha=0.5) + 
  geom_smooth(data=even_points_transformed, aes(x=dist_up_creek, y=dem, group=creek, col=creek), size=0.5) + 
  geom_vline(data=geologic_boundaries, aes(xintercept=(position_rattlesnake+position_mission)/2), col="gray", linetype="dashed") + 
  geom_vline(data=mission_tunnel, aes(xintercept=position), col="magenta", alpha=0.5) +
  ggtitle("Creek Elevation by Longitudinal Position") + 
  scale_color_manual(name="",values = c("mission"="orange","rattlesnake"="green3")) + 
  #scale_linetype_manual(values = c("twodash","dashed")) + 
  xlab("Longitudinal Distance Along Creek (m)") + 
  ylab("Creekbed Elevation (degrees)") + 
  scale_y_continuous(limits=c(-90,900), expand=c(0,0)) + 
  scale_x_continuous(limits=c(0,8000), expand=c(0,0))
ggsave(plot=creek_slope_plot, 
       filename="D:/SERDP/Mission_Creek/output_imagery/mission_rattlesnake_creek_dem.png",
       width=6,height=3,device="png")
creek_slope_plot

#  Slope vs Creek Position
creek_slope_plot <- ggplot() + 
  geom_point(data=even_points_transformed, aes(x=dist_up_creek, y=slope, group=creek, col=creek), alpha=0.5) + 
  geom_smooth(data=even_points_transformed, aes(x=dist_up_creek, y=slope, group=creek, col=creek), size=0.5) + 
  geom_vline(data=geologic_boundaries, aes(xintercept=(position_rattlesnake+position_mission)/2), col="gray", linetype="dashed") + 
  geom_vline(data=mission_tunnel, aes(xintercept=position), col="magenta", alpha=0.5) +
  ggtitle("Creek Slope by Longitudinal Position") + 
  scale_color_manual(name="",values = c("mission"="orange","rattlesnake"="green3")) + 
  #scale_linetype_manual(values = c("twodash","dashed")) + 
  xlab("Longitudinal Distance Along Creek (m)") + 
  ylab("Creekbed Slope (degrees)") + 
  scale_y_continuous(limits=c(-5,55), expand=c(0,0)) + 
  scale_x_continuous(limits=c(0,8000), expand=c(0,0))
ggsave(plot=creek_slope_plot, 
       filename="D:/SERDP/Mission_Creek/output_imagery/mission_rattlesnake_creek_slopes.png",
       width=6,height=3,device="png")
creek_slope_plot

# CHM vs Creek Position
creek_chm_plot <- ggplot(even_points_transformed) + 
  geom_point(aes(x=dist_up_creek, y=chm, group=creek, col=creek), alpha=0.5) + 
  geom_smooth(aes(x=dist_up_creek, y=chm, group=creek, col=creek), size=0.5) + 
  ggtitle("Vegetation Height by Longitudinal Position") + 
  scale_color_manual(name="",values = c("mission"="green3","rattlesnake"="orange")) + 
  #scale_linetype_manual(values = c("twodash","dashed")) + 
  xlab("Longitudinal Distance Along Creek (m)") + 
  ylab("Vegetation Height Over Channel (m)") + 
  scale_y_continuous(limits=c(0,30))
ggsave(plot=creek_chm_plot, 
       filename="D:/SERDP/Mission_Creek/output_imagery/mission_rattlesnake_chm_slopes.png",
       width=6,height=3,device="png")
ggplot(even_points_transformed) + 
  geom_point(aes(x=dist_up_creek, y=slope, group=creek, col=creek)) + 
  geom_smooth(aes(x=dist_up_creek, y=chm, group=creek, col=creek)) + 
  ggtitle("Vegetation Height by Longitudinal Position") + 
  scale_y_continuous(limits=c(0,30)) + 
  xlab("Longitudinal Distance Along Creek (m)") + 
  ylab("Vegetation Height (m)")


# Surface Water vs Creek Position
creek_winter_water_plot <- ggplot(even_points_transformed) + 
  geom_point(aes(x=dist_up_creek, y=dem, group=as.character(!winter_dry), col=as.character(!winter_dry)), alpha=0.5) + 
  ggtitle("Surface Water Distribution in Winter by Longitudinal Position") + 
  scale_color_manual(name="",values = c("FALSE"="red","TRUE"="blue")) + 
  #scale_linetype_manual(values = c("twodash","dashed")) + 
  facet_wrap(~str_to_title(creek), ncol=1) + 
  xlab("Longitudinal Distance Along Creek (m)") + 
  ylab("Creekbed Elevation (m)") + 
  scale_y_continuous(limits=c(-90,900), expand=c(0,0)) + 
  scale_x_continuous(limits=c(0,8000), expand=c(0,0)) + 
  geom_vline(data=geologic_boundaries, aes(xintercept=(position_rattlesnake+position_mission)/2), col="gray", linetype="dashed") + 
  geom_vline(data=mission_tunnel, aes(xintercept=position), col="magenta", alpha=0.5)
ggsave(plot=creek_winter_water_plot, 
       filename="D:/SERDP/Mission_Creek/output_imagery/mission_rattlesnake_winter_water.png",
       width=6,height=6,device="png")
creek_winter_water_plot


# Surface Water vs Creek Position
creek_summer_water_plot <- ggplot(even_points_transformed) + 
  geom_point(aes(x=dist_up_creek, y=dem, group=as.character(!summer_dry), col=as.character(!summer_dry)), alpha=0.5) + 
  ggtitle("Surface Water Distribution in Summer by Longitudinal Position") + 
  scale_color_manual(name="",values = c("FALSE"="red","TRUE"="blue")) + 
  #scale_linetype_manual(values = c("twodash","dashed")) + 
  facet_wrap(~str_to_title(creek), ncol=1) + 
  xlab("Longitudinal Distance Along Creek (m)") + 
  ylab("Creekbed Elevation (m)") + 
  scale_y_continuous(limits=c(-90,900), expand=c(0,0)) + 
  scale_x_continuous(limits=c(0,8000), expand=c(0,0)) + 
  geom_vline(data=geologic_boundaries, aes(xintercept=(position_rattlesnake+position_mission)/2), col="gray", linetype="dashed") + 
  geom_vline(data=mission_tunnel, aes(xintercept=position), col="magenta", alpha=0.5)
ggsave(plot=creek_summer_water_plot, 
       filename="D:/SERDP/Mission_Creek/output_imagery/mission_rattlesnake_summer_water.png",
       width=6,height=6,device="png")
creek_summer_water_plot



# replace NA in surface water values for field_points with 0
field_points[which(is.na(field_points$surface_flow)),]$surface_flow <- FALSE
field_points[which(is.na(field_points$pooling_water)),]$pooling_water <- FALSE
field_points[which(is.na(field_points$seeps_or_springs)),]$seeps_or_springs <- FALSE
# Add "wet" variable for presence of any surface water
field_points$wet <- (field_points$surface_flow + field_points$pooling_water + field_points$seeps_or_springs) > 0 
field_points <- st_transform(field_points, crs(datastack_aq))

# Subset field points to specific regions, for comparison to raster data:
field_points_aq <-  field_points %>% filter(creek_name == "BR")
field_points_copr <-  field_points %>% filter(creek_name %in% c("NCOS","EW"))
field_points_sr <-  field_points %>% filter(creek_name == "SR")
field_points_ac_myc <-  field_points %>% filter(creek_name %in% c("AC","MYC"))
field_points_sym <-  field_points %>% filter(!(creek_name %in% c("BR","EW","NCOS","SR","AC","MYC")))



# Sample rasters in radius around points
#   LiDAR
sample_radius <- 30
aq_samples <- raster::extract(readAll(stack(lidar_aq, crop(datastack_aq,lidar_aq))), field_points_aq, method='simple', buffer=sample_radius, cellnumbers=TRUE) #, fun=mean, na.rm=TRUE)
copr_samples <- raster::extract(readAll(stack(lidar_copr, crop(datastack_copr,lidar_copr))), field_points_copr, method='simple', buffer=sample_radius, cellnumbers=TRUE)
sr_samples <- raster::extract(readAll(stack(lidar_sr, crop(datastack_sr,lidar_sr))), field_points_sr, method='simple', buffer=sample_radius, cellnumbers=TRUE)
ac_myc_samples <- raster::extract(readAll(stack(lidar_ac_myc, crop(datastack_sym,lidar_ac_myc))), field_points_ac_myc, method='simple', buffer=sample_radius, cellnumbers=TRUE)
sym_samples <- raster::extract(readAll(stack(lidar_sym, crop(datastack_sym,lidar_sym))), field_points_sym, method='simple', buffer=sample_radius, cellnumbers=TRUE)
mc_rc_dense <- raster::extract(readAll(stack(lidar_sym, crop(datastack_sym,lidar_sym))), even_points_transformed, method='simple', buffer=sample_radius, cellnumbers=TRUE)
densiometer_samples <- raster::extract(readAll(stack(lidar_sym, crop(datastack_sym,lidar_sym))), densiometer_measurements, method='simple', cellnumbers=TRUE, na.rm=TRUE)
veg_polygon_samples <- raster::extract(readAll(stack(lidar_sym, crop(datastack_sym,lidar_sym))), st_transform(veg_polygons, st_crs(lidar_sym)), method='simple', cellnumbers=TRUE) #, fun=mean, na.rm=TRUE)

# add data to densiometer and veg samples
rowSDs <- function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}
# add data to densiometer and veg samples
densiometer_df <- cbind(densiometer_measurements, densiometer_samples) %>%
  mutate(deciduousness = (sym_fr_datastack_2021.4+sym_fr_datastack_2021.5-sym_fr_datastack_2021.1-sym_fr_datastack_2021.2),
         max_greenness = max(sym_fr_datastack_2021.1, sym_fr_datastack_2021.2, sym_fr_datastack_2021.3, sym_fr_datastack_2021.4, sym_fr_datastack_2021.5, sym_fr_datastack_2021.6, sym_fr_datastack_2021.7, sym_fr_datastack_2021.8, sym_fr_datastack_2021.8, sym_fr_datastack_2021.9, sym_fr_datastack_2021.10, sym_fr_datastack_2021.11, sym_fr_datastack_2021.12),
         min_greenness = max(sym_fr_datastack_2021.1, sym_fr_datastack_2021.2, sym_fr_datastack_2021.3, sym_fr_datastack_2021.4, sym_fr_datastack_2021.5, sym_fr_datastack_2021.6, sym_fr_datastack_2021.7, sym_fr_datastack_2021.8, sym_fr_datastack_2021.8, sym_fr_datastack_2021.9, sym_fr_datastack_2021.10, sym_fr_datastack_2021.11, sym_fr_datastack_2021.12),
         mean_greenness = rowMeans(cbind(sym_fr_datastack_2021.1 + sym_fr_datastack_2021.2+ sym_fr_datastack_2021.3, sym_fr_datastack_2021.4, sym_fr_datastack_2021.5, sym_fr_datastack_2021.6, sym_fr_datastack_2021.7, sym_fr_datastack_2021.8, sym_fr_datastack_2021.8, sym_fr_datastack_2021.9, sym_fr_datastack_2021.10, sym_fr_datastack_2021.11, sym_fr_datastack_2021.12)),
         sd_greenness = rowSDs(cbind(sym_fr_datastack_2021.1, sym_fr_datastack_2021.2, sym_fr_datastack_2021.3, sym_fr_datastack_2021.4, sym_fr_datastack_2021.5, sym_fr_datastack_2021.6, sym_fr_datastack_2021.7, sym_fr_datastack_2021.8, sym_fr_datastack_2021.8, sym_fr_datastack_2021.9, sym_fr_datastack_2021.10, sym_fr_datastack_2021.11, sym_fr_datastack_2021.12)),
         tall_vegetation = (vegetation_histogram_mosaic.5 + vegetation_histogram_mosaic.6 +
                              vegetation_histogram_mosaic.7 + vegetation_histogram_mosaic.8 +
                              vegetation_histogram_mosaic.9 + vegetation_histogram_mosaic.10 +
                              vegetation_histogram_mosaic.11 + vegetation_histogram_mosaic.12 +
                              vegetation_histogram_mosaic.13 + vegetation_histogram_mosaic.14 +
                              vegetation_histogram_mosaic.15 + vegetation_histogram_mosaic.16 +
                              vegetation_histogram_mosaic.17 + vegetation_histogram_mosaic.18 +
                              vegetation_histogram_mosaic.19 + vegetation_histogram_mosaic.20 +
                              vegetation_histogram_mosaic.21 + vegetation_histogram_mosaic.22 +
                              vegetation_histogram_mosaic.23 + vegetation_histogram_mosaic.24 +
                              vegetation_histogram_mosaic.25 + vegetation_histogram_mosaic.26 +
                              vegetation_histogram_mosaic.27 + vegetation_histogram_mosaic.28 +
                              vegetation_histogram_mosaic.29 + vegetation_histogram_mosaic.30 +
                              vegetation_histogram_mosaic.31 + vegetation_histogram_mosaic.32))

overall_cover_fit <- (lm(data=densiometer_df, overall_cover ~ mean_greenness+sd_greenness+deciduousness+tall_vegetation))
summary(overall_cover_fit)
plot(densiometer_df$overall_cover, overall_cover_fit$fitted.values)

decid_cover_fit <- (lm(data=densiometer_df %>% filter(FID > 65), decid_cover ~ mean_greenness+sd_greenness+deciduousness+tall_vegetation))
summary(decid_cover_fit)
plot(densiometer_df[densiometer_df$FID > 65,]$decid_cover, decid_cover_fit$fitted.values)


veg_polygon_df <- cbind(veg_polygons, veg_polygon_samples)

#     Combine creeks into one list
vegetation_samples <- c(aq_samples, copr_samples, sr_samples, ac_myc_samples, sym_samples)
sample_site_names <- as.character(rbind(field_points_aq,
                                        field_points_copr,
                                        field_points_sr,
                                        field_points_ac_myc,
                                        field_points_sym)$site_name)
sample_creek_names <- as.character(rbind(field_points_aq,
                                         field_points_copr,
                                         field_points_sr,
                                         field_points_ac_myc,
                                         field_points_sym)$creek_name)
sample_water_presence <- as.character(rbind(field_points_aq,
                                            field_points_copr,
                                            field_points_sr,
                                            field_points_ac_myc,
                                            field_points_sym)$wet)
# Add site name to vegetation samples around field points
addSiteNameToSamples <- function(index)
{
  output_samples <- as.data.frame(vegetation_samples[[index]])
  names(output_samples) <- c("cell_number","lidar_histogram_lower_bound","lidar_histogram_upper_bound",
                             paste("bin_",1:30,sep=""),
                             paste("month_",1:12,sep=""),
                             "summer_ndvi","winter_ndvi","stdev_ndvi","max_ndvi",
                             "hod_1","hod_2","flowaccum","dem","slope","precip","tmax","tmin","vpd","latitude",
                             "pheno_change_angle","blue","green","red","NIR","SWIR_1","SWIR_2")
  output_samples$site_name <- sample_site_names[[index]]
  output_samples$creek_name <- sample_creek_names[[index]]
  output_samples$wet <- sample_water_presence[[index]]
  return(output_samples)
}
# Add site name to vegetation samples around evenly spaced Mission/Rattlesnake points
addCreekNameToMCRC <- function(index)
{
  output_samples <- as.data.frame(mc_rc_dense[[index]])
  names(output_samples) <- c("cell_number","lidar_histogram_lower_bound","lidar_histogram_upper_bound",
                             paste("bin_",1:30,sep=""),
                             paste("month_",1:12,sep=""),
                             "summer_ndvi","winter_ndvi","stdev_ndvi","max_ndvi",
                             "hod_1","hod_2","flowaccum","dem","slope","precip","tmax","tmin","vpd","latitude",
                             "pheno_change_angle","blue","green","red","NIR","SWIR_1","SWIR_2")
  output_samples$creek_name <- even_points_transformed[index,]$creek
  output_samples$summer_dry <- even_points_transformed[index,]$summer_dry
  output_samples$winter_dry <- even_points_transformed[index,]$winter_dry
  output_samples$winter_water_dist <- even_points_transformed[index,]$winter_water_dist
  output_samples$summer_water_dist <- even_points_transformed[index,]$summer_water_dist
  output_samples$any_pool_dist <- even_points_transformed[index,]$any_pool_dist
  output_samples$creekbed_elevation <- even_points_transformed[index,]$dem
  output_samples$creekbed_slope <- even_points_transformed[index,]$slope
  output_samples$creekbed_chm <- even_points_transformed[index,]$chm
  output_samples$dist_up_creek <- even_points_transformed[index,]$dist_up_creek
  output_samples$formation <- even_points_transformed[index,]$formation
  return(output_samples)
}

vegetation_df <- bind_rows(lapply(1:length(vegetation_samples), addSiteNameToSamples))
mc_rc_df <- bind_rows(lapply(1:length(mc_rc_dense), addCreekNameToMCRC))

modal_value <- function(vector_of_data)
{
  data_df <- data.frame(value = vector_of_data)
  sorted_tally <- data_df %>% 
    group_by(value) %>%
    tally() %>%
    arrange(-n)
  return(sorted_tally[1,]$value)
}

sample_topography <- vegetation_df %>%
  group_by(site_name) %>%
  summarize(elevation = mean(dem),
            slope_min = min(slope),
            slope_max = max(slope),
            slope_median = median(slope), 
            flow_accum = max(flowaccum),
            precip = max(precip),
            tmin = min(tmin),
            tmax = max(tmax),
            vpd = min(vpd),
            creek_name = creek_name[[1]],
            decid_frac = sum((month_6+month_7-month_1-month_2)>0)/n())
sample_topography_mc_rc <- mc_rc_df %>%
  group_by(dist_up_creek, creek_name) %>%
  summarize(elevation = mean(dem),
            formation = modal_value(formation),
            creekbed_elevation = mean(creekbed_elevation),
            creekbed_slope = mean(creekbed_slope),
            creekbed_chm = mean(creekbed_chm),
            summer_dry = mean(summer_dry),
            winter_dry = mean(winter_dry),
            winter_water_dist = mean(winter_water_dist),
            summer_water_dist = mean(summer_water_dist),
            any_pool_dist = mean(any_pool_dist),
            dist_up_creek = mean(dist_up_creek),
            slope_min = min(slope),
            slope_max = max(slope),
            slope_median = median(slope), 
            flow_accum = max(flowaccum),
            precip = max(precip),
            tmin = min(tmin),
            tmax = max(tmax),
            vpd = min(vpd),
            creek_name = creek_name[[1]],
            decid_frac = sum((month_6+month_7-month_1-month_2)>0)/n(),
            max_greenness = max(month_1,month_2,month_3,month_4,month_5,month_6,month_7,month_8,month_9,month_10,month_11,month_12),
            min_greenness = min(month_1,month_2,month_3,month_4,month_5,month_6,month_7,month_8,month_9,month_10,month_11,month_12),
            median_greenness = median(month_1,month_2,month_3,month_4,month_5,month_6,month_7,month_8,month_9,month_10,month_11,month_12),
            mean_return_height = mean(bin_1*1+bin_2*2+bin_3*3+bin_4*4+bin_5*5+bin_6*6+bin_7*7+bin_8*8+bin_9*9+bin_10*10+bin_11*11+bin_12*12+bin_13*13+bin_14*14+bin_15*15+bin_16*16+bin_17*17+bin_18*18+bin_19*19+bin_20*20+bin_21*21+bin_22*22+bin_23*23+bin_24*24+bin_25*25+bin_26*26+bin_27*27+bin_28*28+bin_29*29+bin_30*30),
            max_return_height = max((bin_1>0)*1,(bin_2>0)*2,(bin_3>0)*3,(bin_4>0)*4,(bin_5>0)*5,(bin_6>0)*6,(bin_7>0)*7,(bin_8>0)*8,(bin_9>0)*9,(bin_10>0)*10,(bin_11>0)*11,(bin_12>0)*12,(bin_13>0)*13,(bin_14>0)*14,(bin_15>0)*15,(bin_16>0)*16,(bin_17>0)*17,(bin_18>0)*18,(bin_19>0)*19,(bin_20>0)*20,(bin_21>0)*21,(bin_22>0)*22,(bin_23>0)*23,(bin_24>0)*24,(bin_25>0)*25,(bin_26>0)*26,(bin_27>0)*27,(bin_28>0)*28,(bin_29>0)*29,(bin_30>0)*30),
            min_return_height = min((bin_1>0)*1,(bin_2>0)*2,(bin_3>0)*3,(bin_4>0)*4,(bin_5>0)*5,(bin_6>0)*6,(bin_7>0)*7,(bin_8>0)*8,(bin_9>0)*9,(bin_10>0)*10,(bin_11>0)*11,(bin_12>0)*12,(bin_13>0)*13,(bin_14>0)*14,(bin_15>0)*15,(bin_16>0)*16,(bin_17>0)*17,(bin_18>0)*18,(bin_19>0)*19,(bin_20>0)*20,(bin_21>0)*21,(bin_22>0)*22,(bin_23>0)*23,(bin_24>0)*24,(bin_25>0)*25,(bin_26>0)*26,(bin_27>0)*27,(bin_28>0)*28,(bin_29>0)*29,(bin_30>0)*30))

vegetation_structure_long <- vegetation_df %>% 
  select(c(paste("bin_",1:30,sep=""),"site_name","creek_name","wet")) %>%
  pivot_longer(1:30,names_to="bin_string",values_to="fraction") %>%
  mutate(bin = as.numeric(substring(bin_string, 5,10)))
vegetation_structure_long <- merge(vegetation_structure_long,(sample_topography %>% select(-creek_name)),by="site_name")
vegetation_phenology_long <- vegetation_df %>% 
  select(c(paste("month_",1:12,sep=""),"site_name","creek_name","wet")) %>%
  pivot_longer(1:12,names_to="month_string",values_to="NDVI") %>%
  mutate(month = as.numeric(substring(month_string, 7,10)))
vegetation_phenology_long <- merge(vegetation_phenology_long,(sample_topography %>% select(-creek_name)),by="site_name")

vegetation_structure_long_mcrc <- mc_rc_df %>% 
  select(c(paste("bin_",1:30,sep=""),"creek_name","dem","winter_dry","summer_dry","winter_water_dist","summer_water_dist","any_pool_dist","creekbed_elevation","creekbed_slope","creekbed_chm","dist_up_creek")) %>%
  mutate(elevation=dem) %>%
  pivot_longer(1:30,names_to="bin_string",values_to="fraction") %>%
  mutate(bin = as.numeric(substring(bin_string, 5,10)))
vegetation_structure_long_mcrc <- merge(vegetation_structure_long_mcrc,(sample_topography_mc_rc %>% select(-creek_name)),by="elevation")
vegetation_phenology_long_mcrc <- mc_rc_df %>% 
  select(c(paste("month_",1:12,sep=""),"creek_name","dem","winter_dry","summer_dry","winter_water_dist","summer_water_dist","any_pool_dist","creekbed_elevation","creekbed_slope","creekbed_chm","dist_up_creek")) %>%
  mutate(elevation=dem) %>%
  pivot_longer(1:12,names_to="month_string",values_to="NDVI") %>%
  mutate(month = as.numeric(substring(month_string, 7,10)))
vegetation_phenology_long_mcrc <- merge(vegetation_phenology_long_mcrc,(sample_topography_mc_rc %>% select(-creek_name)),by="elevation")



# Visualize structure across Rattlesnake and Mission Creek sites - based on averages 
mc_rc_structure_centers <- ggplot(vegetation_structure_long_mcrc %>% drop_na(formation, creek_name) %>% mutate(creek_name = str_to_title(creek_name))) + 
  geom_smooth(aes(y=fraction,x=bin, group=creek_name, col=creek_name), se=FALSE) + 
  facet_wrap(~formation, nrow=1) + 
  scale_y_continuous(limits=c(.01,0.25), expand=c(0,0)) + 
  scale_x_continuous(limits=c(0,25), expand=c(0,0)) + 
  scale_fill_manual(values=c(colorRampPalette(paletteer::paletteer_dynamic("cartography::green.pal", 
                                                                           n = 20, 
                                                                           direction = 1))(10)), 
                    aesthetics = c("fill", "color")) + 
  ylab("Fraction of Returns") + 
  xlab("Vegetation Height") + 
  ggtitle("Vegetation Height Distribution by Creek and Geologic Formation") + 
  coord_flip() + 
  scale_color_manual(name="",values = c("Mission"="orange","Rattlesnake"="green3")) 
mc_rc_structure_centers
ggsave("D:/SERDP/Mission_Creek/output_imagery/mc_rc_lidar_structure_mean.png", mc_rc_structure_centers, width=8, height=3)
# Visualize structure across Rattlesnake and Mission Creek sites - based on averages 
mc_rc_structure_density <- ggplot(vegetation_structure_long_mcrc %>% drop_na(formation, creek_name) %>% mutate(creek_name = str_to_title(creek_name))) + 
  geom_density_2d_filled(aes(y=fraction, x=bin), contour_var = "ndensity") + 
  facet_wrap(~creek_name+formation, nrow=2) + 
  scale_y_continuous(limits=c(.01,0.3), expand=c(0,0)) + 
  scale_x_continuous(limits=c(0,25), expand=c(0,0)) + 
  scale_fill_manual(values=c(colorRampPalette(paletteer::paletteer_dynamic("cartography::green.pal", 
                                                                           n = 20, 
                                                                           direction = 1))(10)), 
                    aesthetics = c("fill", "color")) + 
  ylab("Fraction of Returns") + 
  xlab("Vegetation Height") + 
  ggtitle("Vegetation Height Distribution by Creek and Geologic Formation") + 
  coord_flip()
mc_rc_structure_density
ggsave("D:/SERDP/Mission_Creek/output_imagery/mc_rc_lidar_structure_density.png", mc_rc_structure_density, width=12, height=4)



# Visualize structure across Rattlesnake and Mission Creek sites
mc_rc_phenology_centers <- ggplot(vegetation_phenology_long_mcrc %>% drop_na(creek_name, formation)) + 
  geom_smooth(aes(x=month, y=NDVI, col=str_to_title(creek_name), group=str_to_title(creek_name)), se=FALSE) + 
  facet_wrap(~formation, nrow=1) + 
  scale_y_continuous(limits=c(0.5,1)) +
  scale_x_continuous(limits=c(1,12), expand=c(0,0)) + 
  ylab("Normalized Difference Vegetation Index") + 
  xlab("Vegetation Height") + 
  ggtitle("Riparian Woodland Leaf Phenology by Creek and Geologic Formation") + 
  scale_color_manual(name="",values = c("Mission"="orange","Rattlesnake"="green3")) 
mc_rc_phenology_centers
ggsave("D:/SERDP/Mission_Creek/output_imagery/mc_rc_phenology_mean.png", mc_rc_phenology_centers, width=12, height=4)

ggplot(vegetation_structure_long %>%
         filter(creek_name %in% c("RC","MC"))) + 
  geom_density_2d_filled(aes(y=bin, x=fraction), contour_var = "ndensity") + 
  facet_wrap(~site_name) + 
  scale_x_continuous(limits=c(.01,0.3)) + 
  scale_fill_manual(values=c(colorRampPalette(paletteer::paletteer_dynamic("cartography::green.pal", 
                                                                           n = 20, 
                                                                           direction = 1))(10)), 
                    aesthetics = c("fill", "color"))


# Visualize structure across Rattlesnake and Mission Creek sites
mc_rc_phenology_centers <- ggplot(vegetation_phenology_long_mcrc %>% drop_na(creek_name, formation)) + 
  geom_smooth(aes(x=month, y=NDVI, col=str_to_title(creek_name), group=str_to_title(creek_name)), se=FALSE) + 
  facet_wrap(~formation, nrow=1) + 
  scale_y_continuous(limits=c(0.5,1)) +
  scale_x_continuous(limits=c(1,12), expand=c(0,0)) + 
  ylab("Normalized Difference Vegetation Index") + 
  xlab("Vegetation Height") + 
  ggtitle("Riparian Woodland Leaf Phenology by Creek and Geologic Formation") + 
  scale_color_manual(name="",values = c("Mission"="orange","Rattlesnake"="green3")) 
mc_rc_phenology_centers
ggsave("D:/SERDP/Mission_Creek/output_imagery/mc_rc_phenology_mean.png", mc_rc_phenology_centers, width=12, height=4)
# Density plots 
mc_rc_phenology_density <- ggplot(vegetation_phenology_long_mcrc %>% drop_na(creek_name, formation)) + 
  geom_density_2d_filled(aes(x=month, y=NDVI), contour_var = "ndensity") + 
  geom_smooth(aes(x=month, y=NDVI), se=FALSE, col="red") + 
  facet_wrap(~creek_name+formation, nrow=2) + 
  scale_y_continuous(limits=c(0.5,1)) +
  scale_x_continuous(limits=c(1,12), expand=c(0,0)) + 
  ylab("Normalized Difference Vegetation Index") + 
  xlab("Vegetation Height") + 
  ggtitle("Riparian Woodland Leaf Phenology by Creek and Geologic Formation")
mc_rc_phenology_density
ggsave("D:/SERDP/Mission_Creek/output_imagery/mc_rc_phenology_density.png", mc_rc_phenology_density, width=12, height=4)


# Unsupervised clustering on LiDAR + phenology, binning ALL sampled cells together
set.seed(3)
num_clusters <- 10
veg_clusters <- stats::kmeans(na.omit((vegetation_df)[,4:45]), centers=num_clusters, nstart=10)
veg_cluster_df <- as.data.frame(veg_clusters$centers) %>%
  mutate(cluster = 1:num_clusters) %>%
  pivot_longer(1:42, names_to="bin", values_to="fraction") %>%
  mutate(bin_num = as.numeric(unlist(lapply(str_split(bin,"_"),function(x){return(x[[2]])}))) + 
           as.numeric(unlist(lapply(str_split(bin,"_"),function(x){return(x[[1]])}))=="month")*30)
ggplot(veg_cluster_df) + 
  geom_line(aes(x=bin_num, y=fraction, col=as.factor(cluster), group=cluster))
ggplot(veg_cluster_df) + 
  geom_line(aes(x=bin_num, y=fraction*(1+(bin_num<=30)*2), col=bin_num>30)) + 
  facet_wrap(~cluster)

getHeightQuantile <- function(q)
{
  veg_collapsed <- vegetation_df %>% 
    group_by(site_name) %>%
    summarize(bin_1 = mean(bin_1),
              bin_2 = mean(bin_2),
              bin_3 = mean(bin_3),
              bin_4 = mean(bin_4),
              bin_5 = mean(bin_5),
              bin_6 = mean(bin_6),
              bin_7 = mean(bin_7),
              bin_8 = mean(bin_8),
              bin_9 = mean(bin_9),
              bin_10 = mean(bin_10),
              bin_11 = mean(bin_11),
              bin_12 = mean(bin_12),
              bin_13 = mean(bin_13),
              bin_14 = mean(bin_14),
              bin_15 = mean(bin_15),
              bin_16 = mean(bin_16),
              bin_17 = mean(bin_17),
              bin_18 = mean(bin_18),
              bin_19 = mean(bin_19),
              bin_20 = mean(bin_20),
              bin_21 = mean(bin_21),
              bin_22 = mean(bin_22),
              bin_23 = mean(bin_23),
              bin_24 = mean(bin_24),
              bin_25 = mean(bin_25),
              bin_26 = mean(bin_26),
              bin_27 = mean(bin_27),
              bin_28 = mean(bin_28),
              bin_29 = mean(bin_29),
              bin_30 = mean(bin_30))
  return(as.numeric(unlist(lapply(1:nrow(veg_collapsed), 
                           FUN=function(y){
                             quantile(unlist(lapply(1:30, 
                                                    FUN=function(x){ 
                                                      return(rep(x, round((as.matrix(veg_collapsed[y,x+1]))*100)))
                                                      })), 
                                      q)}))))
}

vegetation_df$veg_cluster <- veg_clusters$cluster
veg_df_collapsed <- vegetation_df %>% 
  mutate(max_height = apply((as.matrix(vegetation_df[,4:33])>0) %*% diag((1:30)), 1, max),
         mean_height = apply((as.matrix(vegetation_df[,4:33])) %*% diag((1:30)), 1, sum)) %>%
  group_by(site_name) %>% 
  summarize(cluster_1 = mean(veg_cluster == 1),
            cluster_2 = mean(veg_cluster == 2),
            cluster_3 = mean(veg_cluster == 3), 
            cluster_4 = mean(veg_cluster == 4), 
            cluster_5 = mean(veg_cluster == 5),
            cluster_6 = mean(veg_cluster == 6),
            cluster_7 = mean(veg_cluster == 7),
            cluster_8 = mean(veg_cluster == 8),
            cluster_9 = mean(veg_cluster == 9),
            cluster_10 = mean(veg_cluster == 10),
            low_veg = mean(bin_1+bin_2),
            mid_veg = mean(bin_3+bin_4+bin_5),
            tall_veg = mean(bin_6 + bin_7 + bin_8 +
               bin_9 + bin_10 +
               bin_11 + bin_12 +
               bin_13 + bin_14 +
               bin_15 + bin_16 +
               bin_17 + bin_18 +
               bin_19 + bin_20 +
               bin_21 + bin_22 +
               bin_23 + bin_24 +
               bin_25 + bin_26 +
               bin_27 + bin_28 +
               bin_29 + bin_30 ),
            max_height = max(max_height),
            mean_height = mean(mean_height),
            mean_max_height = mean(max_height),
            max_over_mean_height = max_height / mean_height,
            decid = sum(summer_ndvi > winter_ndvi)/n(),
            summer_greenness = sum(summer_ndvi*(summer_ndvi>0.4), na.rm=TRUE)/sum((summer_ndvi>0.4), na.rm=TRUE),
            dense_veg_fraction = sum(summer_ndvi > 0.4)/n(),
            std_max_greenness = sd(max_ndvi),
            mean_std_greenness = mean(stdev_ndvi),
            elevation = mean(dem),
            std_elevation = sd(dem)) 
veg_df_collapsed$height_pct_95 <- getHeightQuantile(0.95)/veg_df_collapsed$max_height
veg_df_collapsed$height_pct_90 <- getHeightQuantile(0.9)/veg_df_collapsed$max_height
veg_df_collapsed$height_pct_80 <- getHeightQuantile(0.8)/veg_df_collapsed$max_height
veg_df_collapsed$height_pct_70 <- getHeightQuantile(0.7)/veg_df_collapsed$max_height
veg_df_collapsed$height_pct_60 <- getHeightQuantile(0.6)/veg_df_collapsed$max_height
veg_df_collapsed$height_pct_50 <- getHeightQuantile(0.5)/veg_df_collapsed$max_height
veg_df_collapsed$height_pct_40 <- getHeightQuantile(0.4)/veg_df_collapsed$max_height
veg_df_collapsed$height_pct_30 <- getHeightQuantile(0.3)/veg_df_collapsed$max_height
veg_df_collapsed$height_pct_20 <- getHeightQuantile(0.2)/veg_df_collapsed$max_height
veg_df_collapsed$height_pct_10 <- getHeightQuantile(0.1)/veg_df_collapsed$max_height
veg_df_collapsed$summer_greenness <- replace(veg_df_collapsed$summer_greenness, is.na(veg_df_collapsed$summer_greenness), 0)
sites_with_audio <- unique(all_summaries$site_name)
veg_with_audio <- veg_df_collapsed %>% 
  filter(site_name %in% sites_with_audio)
# YEWA Presence Model
vector_mode <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}
creek_from_site <- all_summaries %>% 
  group_by(site_name) %>%
  summarize(creek = vector_mode(creek),
            site_name_x = vector_mode(site_name))
all_summaries_collapsed <- all_summaries %>%
  group_by(common_name, site_name) %>% 
  summarize(best_confidence = max(best_confidence),
            count = sum(count),
            wet = max(wet),
            creek = vector_mode(creek))

# Normalize predictor variables: 
veg_with_audio <- veg_with_audio %>% 
  mutate(max_height = (max_height-min(max_height))/(max(max_height)-min(max_height)),
         mean_height = (mean_height-min(mean_height))/(max(mean_height)-min(mean_height)),
         mean_max_height = (mean_max_height-min(mean_max_height))/(max(mean_max_height)-min(mean_max_height)),
         max_over_mean_height = (max_over_mean_height-min(max_over_mean_height))/(max(max_over_mean_height)-min(max_over_mean_height)),
         decid = (decid-min(decid))/(max(decid)-min(decid)),
         summer_greenness = (summer_greenness-min(summer_greenness))/(max(summer_greenness)-min(summer_greenness)),
         dense_veg_fraction = (dense_veg_fraction-min(dense_veg_fraction))/(max(dense_veg_fraction)-min(dense_veg_fraction)),
         std_max_greenness = (std_max_greenness-min(std_max_greenness))/(max(std_max_greenness)-min(std_max_greenness)),
         mean_std_greenness = (mean_std_greenness-min(mean_std_greenness))/(max(mean_std_greenness)-min(mean_std_greenness)),
         elevation = (elevation-min(elevation))/(max(elevation)-min(elevation)),
         std_elevation = (std_elevation-min(std_elevation))/(max(std_elevation)-min(std_elevation)),
         low_veg = (low_veg-min(low_veg))/(max(low_veg)-min(low_veg)),
         mid_veg = (mid_veg-min(mid_veg))/(max(mid_veg)-min(mid_veg)),
         tall_veg = (tall_veg-min(tall_veg))/(max(tall_veg)-min(tall_veg)))

veg_audio_summaries_yewa <- merge(veg_with_audio, all_summaries_collapsed %>% filter(common_name == "Yellow Warbler"), by="site_name", all=TRUE) %>%
  mutate(present = !is.na(count)) %>%
  mutate(creek_name = creek_from_site[creek_from_site$site_name_x==site_name,]$creek)
veg_audio_summaries_yewa <- replace(veg_audio_summaries_yewa, is.na(veg_audio_summaries_yewa), 0)
veg_audio_summaries_yewa <- veg_audio_summaries_yewa %>%
  mutate(present_restricted = present*(count>50)*(best_confidence>0.85))

sum(veg_audio_summaries_yewa$present_restricted)
veg_audio_summaries_yewa %>% filter(present_restricted == 1) %>% select(creek, site_name, decid, tall_veg, best_confidence, count) %>% arrange(best_confidence)

logistic_model_yewa <- glm(present_restricted ~ decid+summer_greenness+std_max_greenness+mean_std_greenness+tall_veg+low_veg+wet+elevation+std_elevation+max_height+max_over_mean_height, family = binomial(), data=veg_audio_summaries_yewa)
summary(logistic_model_yewa)
logistic_model_yewa <- glm(present_restricted ~ decid+dense_veg_fraction+wet, family = binomial(), data=veg_audio_summaries_yewa)
summary(logistic_model_yewa)

veg_audio_yewa_sf <- merge(field_data, veg_audio_summaries_yewa, by="site_name") %>% st_as_sf(sf_column_name = "geometry")
st_write(obj=veg_audio_yewa_sf, dsn=here::here("vegetation_and_audio_yewa.gpkg"), delete_dsn = TRUE)




# WAVI Presence Model
veg_audio_summaries_wavi <- merge(veg_with_audio, all_summaries_collapsed %>% filter(common_name == "Warbling Vireo"), by="site_name", all=TRUE) %>%
  mutate(present = !is.na(count)) %>%
  mutate(creek_name = creek_from_site[creek_from_site$site_name_x==site_name,]$creek)
veg_audio_summaries_wavi <- replace(veg_audio_summaries_wavi, is.na(veg_audio_summaries_wavi), 0)
veg_audio_summaries_wavi <- veg_audio_summaries_wavi %>%
  mutate(present_restricted = present*(count>50)*(best_confidence>0.85))

logistic_model_wavi <- glm(present_restricted ~ decid+summer_greenness+std_max_greenness+mean_std_greenness+tall_veg+low_veg+wet+elevation+std_elevation+max_height+max_over_mean_height, family = binomial(), data=veg_audio_summaries_wavi)
summary(logistic_model_wavi)

logistic_model_wavi <- glm(present_restricted ~ tall_veg+low_veg, family = binomial(), data=veg_audio_summaries_wavi)
summary(logistic_model_wavi)


# WIWA Presence Model
veg_audio_summaries_wiwa <- merge(veg_with_audio, all_summaries_collapsed %>% filter(common_name == "Wilson's Warbler"), by="site_name", all=TRUE) %>%
  mutate(present = !is.na(count)) %>%
  mutate(creek_name = creek_from_site[creek_from_site$site_name_x==site_name,]$creek)
veg_audio_summaries_wiwa <- replace(veg_audio_summaries_wiwa, is.na(veg_audio_summaries_wiwa), 0)
veg_audio_summaries_wiwa <- veg_audio_summaries_wiwa %>%
  mutate(present_restricted = present*(count>50)*(best_confidence>0.85))

logistic_model_wiwa <- glm(present_restricted ~ decid+summer_greenness+std_max_greenness+mean_std_greenness+tall_veg+low_veg+wet+elevation+std_elevation+max_height+max_over_mean_height, family = binomial(), data=veg_audio_summaries_wiwa)
summary(logistic_model_wiwa)

logistic_model_wiwa <- glm(present_restricted ~ wet, family = binomial(), data=veg_audio_summaries_wiwa)
summary(logistic_model_wiwa)



# ACWO Presence Model
veg_audio_summaries_acwo <- merge(veg_with_audio, all_summaries_collapsed %>% filter(common_name == "Acorn Woodpecker"), by="site_name", all=TRUE) %>%
  mutate(present = !is.na(count)) %>%
  mutate(creek_name = creek_from_site[creek_from_site$site_name_x==site_name,]$creek)
veg_audio_summaries_acwo <- replace(veg_audio_summaries_acwo, is.na(veg_audio_summaries_acwo), 0)
veg_audio_summaries_acwo <- veg_audio_summaries_acwo %>%
  mutate(present_restricted = present*(count>50)*(best_confidence>0.85))

logistic_model_acwo <- glm(present_restricted ~ decid+summer_greenness+std_max_greenness+mean_std_greenness+tall_veg+low_veg+wet+elevation+std_elevation+max_height+max_over_mean_height, family = binomial(), data=veg_audio_summaries_acwo)
summary(logistic_model_acwo)

logistic_model_acwo <- glm(present_restricted ~ decid+tall_veg+wet+max_over_mean_height, family = binomial(), data=veg_audio_summaries_acwo)
summary(logistic_model_acwo)


# WBNU Presence Model
veg_audio_summaries_wbnu <- merge(veg_with_audio, all_summaries_collapsed %>% filter(common_name == "White-breasted Nuthatch"), by="site_name", all=TRUE) %>%
  mutate(present = !is.na(count)) %>%
  mutate(creek_name = creek_from_site[creek_from_site$site_name_x==site_name,]$creek)
veg_audio_summaries_wbnu <- replace(veg_audio_summaries_wbnu, is.na(veg_audio_summaries_wbnu), 0)
veg_audio_summaries_wbnu <- veg_audio_summaries_wbnu %>%
  mutate(present_restricted = present*(count>50)*(best_confidence>0.85))

logistic_model_wbnu <- glm(present_restricted ~ decid+summer_greenness+std_max_greenness+mean_std_greenness+tall_veg+low_veg+wet+elevation+std_elevation+max_height+max_over_mean_height, family = binomial(), data=veg_audio_summaries_wbnu)
summary(logistic_model_wbnu)

logistic_model_wbnu <- glm(present_restricted ~ decid+tall_veg+wet+max_over_mean_height, family = binomial(), data=veg_audio_summaries_wbnu)
summary(logistic_model_wbnu)



# PUFI Presence Model
veg_audio_summaries_pufi <- merge(veg_with_audio, all_summaries_collapsed %>% filter(common_name == "Purple Finch"), by="site_name", all=TRUE) %>%
  mutate(present = !is.na(count)) %>%
  mutate(creek_name = creek_from_site[creek_from_site$site_name_x==site_name,]$creek)
veg_audio_summaries_pufi <- replace(veg_audio_summaries_pufi, is.na(veg_audio_summaries_pufi), 0)
veg_audio_summaries_pufi <- veg_audio_summaries_pufi %>%
  mutate(present_restricted = present*(count>50)*(best_confidence>0.85))

logistic_model_pufi <- glm(present_restricted ~ decid+summer_greenness+std_max_greenness+mean_std_greenness+tall_veg+low_veg+wet+elevation+std_elevation+max_height+max_over_mean_height, family = binomial(), data=veg_audio_summaries_pufi)
summary(logistic_model_pufi)


# CAWR Presence Model
veg_audio_summaries_cawr <- merge(veg_with_audio, all_summaries_collapsed %>% filter(common_name == "Canyon Wren"), by="site_name", all=TRUE) %>%
  mutate(present = !is.na(count)) %>%
  mutate(creek_name = creek_from_site[creek_from_site$site_name_x==site_name,]$creek)
veg_audio_summaries_cawr <- replace(veg_audio_summaries_cawr, is.na(veg_audio_summaries_cawr), 0)
veg_audio_summaries_cawr <- veg_audio_summaries_cawr %>%
  mutate(present_restricted = present*(count>50)*(best_confidence>0.85))

logistic_model_cawr <- glm(present_restricted ~ decid+summer_greenness+std_max_greenness+mean_std_greenness+tall_veg+low_veg+wet+elevation+std_elevation+max_height+max_over_mean_height, family = binomial(), data=veg_audio_summaries_cawr)
summary(logistic_model_cawr)

logistic_model_cawr <- glm(present_restricted ~ elevation+wet+max_over_mean_height, family = binomial(), data=veg_audio_summaries_cawr)
summary(logistic_model_cawr)



# PSFL Presence Model
veg_audio_summaries_psfl <- merge(veg_with_audio, all_summaries_collapsed %>% filter(common_name == "Pacific-slope Flycatcher"), by="site_name", all=TRUE) %>%
  mutate(present = !is.na(count)) %>%
  mutate(creek_name = creek_from_site[creek_from_site$site_name_x==site_name,]$creek)
veg_audio_summaries_psfl <- replace(veg_audio_summaries_psfl, is.na(veg_audio_summaries_psfl), 0)
veg_audio_summaries_psfl <- veg_audio_summaries_psfl %>%
  mutate(present_restricted = present*(count>50)*(best_confidence>0.85))

logistic_model_psfl <- glm(present_restricted ~ decid+summer_greenness+std_max_greenness+mean_std_greenness+tall_veg+low_veg+wet+elevation+std_elevation+max_height+max_over_mean_height, family = binomial(), data=veg_audio_summaries_psfl)
summary(logistic_model_psfl)

logistic_model_psfl <- glm(present_restricted ~ summer_greenness+wet+elevation+std_elevation, family = binomial(), data=veg_audio_summaries_psfl)
summary(logistic_model_psfl)


# AMRO Presence Model
veg_audio_summaries_amro <- merge(veg_with_audio, all_summaries_collapsed %>% filter(common_name == "American Robin"), by="site_name", all=TRUE) %>%
  mutate(present = !is.na(count)) %>%
  mutate(creek_name = creek_from_site[creek_from_site$site_name_x==site_name,]$creek)
veg_audio_summaries_amro <- replace(veg_audio_summaries_amro, is.na(veg_audio_summaries_amro), 0)
veg_audio_summaries_amro <- veg_audio_summaries_amro %>%
  mutate(present_restricted = present*(count>50)*(best_confidence>0.85))

logistic_model_amro <- glm(present_restricted ~ decid+summer_greenness+std_max_greenness+mean_std_greenness+tall_veg+low_veg+wet+elevation+std_elevation+max_height+max_over_mean_height, family = binomial(), data=veg_audio_summaries_amro)
summary(logistic_model_amro)

logistic_model_amro <- glm(present_restricted ~ summer_greenness+wet+elevation+std_elevation, family = binomial(), data=veg_audio_summaries_amro)
summary(logistic_model_amro)

# SOSP Presence Model
veg_audio_summaries_sosp <- merge(veg_with_audio, all_summaries_collapsed %>% filter(common_name == "Song Sparrow"), by="site_name", all=TRUE) %>%
  mutate(present = !is.na(count)) %>%
  mutate(creek_name = creek_from_site[creek_from_site$site_name_x==site_name,]$creek)
veg_audio_summaries_sosp <- replace(veg_audio_summaries_sosp, is.na(veg_audio_summaries_sosp), 0)
veg_audio_summaries_sosp <- veg_audio_summaries_sosp %>%
  mutate(present_restricted = present*(count>50)*(best_confidence>0.85))

logistic_model_sosp <- glm(present_restricted ~ decid+summer_greenness+std_max_greenness+mean_std_greenness+tall_veg+low_veg+wet+elevation+std_elevation+max_height+max_over_mean_height, family = binomial(), data=veg_audio_summaries_sosp)
summary(logistic_model_sosp)

logistic_model_sosp <- glm(present_restricted ~ decid+summer_greenness+wet+elevation+max_height, family = binomial(), data=veg_audio_summaries_sosp)
summary(logistic_model_sosp)




# Community stats model
site_summaries <- all_summaries %>% 
  group_by(site_name) %>% 
  summarize(total_counts = sum(count),
            richness = n(),
            wet=max(wet))
site_summaries_veg <- merge(site_summaries, veg_with_audio, by="site_name")
summary(lm(data=site_summaries_veg, richness ~ decid+summer_greenness+wet+elevation+max_height))



# Exploratory plots of vegetative characteristics by elevation in the two creeks 
# Mean vegetation return height (from LiDAR)
veg_mean_height_plot <- ggplot(sample_topography_mc_rc) + 
  geom_smooth(aes(x=dist_up_creek,y=mean_return_height,group=str_to_title(creek_name),col=str_to_title(creek_name)), level=0.95, method="loess") + 
  geom_point(aes(x=dist_up_creek,y=mean_return_height,group=str_to_title(creek_name),col=str_to_title(creek_name)), alpha=0.5) +
  geom_vline(data=geologic_boundaries, aes(xintercept=(position_rattlesnake+position_mission)/2), col="gray", linetype="dashed") + 
  geom_vline(data=mission_tunnel, aes(xintercept=position), col="magenta", alpha=0.5) +
  geom_vline(data=summer_water_boundaries, aes(xintercept=position, col=creek), alpha=0.5, linetype="dashed") +
  ggtitle("Mean Vegetation Height by Longitudinal Position") + 
  scale_y_continuous(limits=c(0,15), expand=c(0,0)) +
  scale_x_continuous(limits=c(0,8000), expand=c(0,0)) +
  xlab("Longitudinal Distance Along Creek (m)") + 
  ylab("Mean Vegetation Height (m)") +
  scale_color_manual(name="",values = c("Mission"="orange","Rattlesnake"="green3")) 
veg_mean_height_plot
ggsave("D:/SERDP/Mission_Creek/output_imagery/mean_veg_height.png", veg_mean_height_plot, width=6, height=3, device="png")
# Max vegetation return height (from LiDAR)
veg_max_height_plot <- ggplot(sample_topography_mc_rc) + 
  geom_smooth(aes(x=dist_up_creek,y=max_return_height,group=str_to_title(creek_name),col=str_to_title(creek_name)), level=0.95, method="loess") + 
  geom_point(aes(x=dist_up_creek,y=max_return_height,group=str_to_title(creek_name),col=str_to_title(creek_name)), alpha=0.5) + 
  geom_vline(data=geologic_boundaries, aes(xintercept=(position_rattlesnake+position_mission)/2), col="gray", linetype="dashed") + 
  geom_vline(data=mission_tunnel, aes(xintercept=position), col="magenta", alpha=0.5) +
  geom_vline(data=summer_water_boundaries, aes(xintercept=position, col=creek), alpha=0.5, linetype="dashed") +
  ggtitle("Maximum Vegetation Height by Longitudinal Position") + 
  scale_y_continuous(limits=c(0,30), expand=c(0,0)) +
  scale_x_continuous(limits=c(0,8000), expand=c(0,0)) +
  xlab("Longitudinal Distance Along Creek (m)") + 
  ylab("Maximum Vegetation Height (m)") +
  scale_color_manual(name="",values = c("Mission"="orange","Rattlesnake"="green3")) 
veg_max_height_plot
ggsave("D:/SERDP/Mission_Creek/output_imagery/max_veg_height.png", veg_max_height_plot, width=6, height=3, device="png")
# Deciduousness (from GEE Phenology)
veg_decid_frac <- ggplot(sample_topography_mc_rc) + 
  geom_smooth(aes(x=dist_up_creek,y=decid_frac,group=str_to_title(creek_name),col=str_to_title(creek_name)), level=0.95, method="loess") + 
  geom_point(aes(x=dist_up_creek,y=decid_frac,group=str_to_title(creek_name),col=str_to_title(creek_name)), alpha=0.5) + 
  geom_vline(data=geologic_boundaries, aes(xintercept=(position_rattlesnake+position_mission)/2), col="gray", linetype="dashed") + 
  geom_vline(data=mission_tunnel, aes(xintercept=position), col="magenta", alpha=0.5) +
  geom_vline(data=summer_water_boundaries, aes(xintercept=position, col=creek), alpha=0.5, linetype="dashed") +
  ggtitle("Deciduous Vegetation Cover by Longitudinal Position") + 
  scale_y_continuous(limits=c(-0.1,1), expand=c(0,0)) +
  scale_x_continuous(limits=c(-0.1,8000), expand=c(0,0)) +
  xlab(" Longitudinal Position (m)") + 
  ylab("Fraction of Winter-deciduous Vegetation") +
  scale_color_manual(name="",values = c("Mission"="orange","Rattlesnake"="green3")) 
veg_decid_frac
ggsave("D:/SERDP/Mission_Creek/output_imagery/veg_decid_frac.png", veg_decid_frac, width=6, height=3, device="png")
# Deciduousness (from GEE Phenology)
# Mean Greenness (from GEE Phenology)
veg_median_greenness <- ggplot(sample_topography_mc_rc) + 
  geom_smooth(aes(x=dist_up_creek,y=median_greenness,group=str_to_title(creek_name),col=str_to_title(creek_name)), level=0.95, method="loess") + 
  geom_point(aes(x=dist_up_creek,y=median_greenness,group=str_to_title(creek_name),col=str_to_title(creek_name)), alpha=0.5) + 
  geom_vline(data=geologic_boundaries, aes(xintercept=(position_rattlesnake+position_mission)/2), col="gray", linetype="dashed") + 
  geom_vline(data=mission_tunnel, aes(xintercept=position), col="magenta", alpha=0.5) +
  geom_vline(data=summer_water_boundaries, aes(xintercept=position, col=creek), alpha=0.5, linetype="dashed") +
  ggtitle("Median Vegetation Greenness by Longitudinal Position") + 
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  scale_x_continuous(limits=c(0,8000), expand=c(0,0)) +
  xlab("Longitudinal Position (m)") + 
  ylab("Median Greenness") +
  scale_color_manual(name="",values = c("Mission"="orange","Rattlesnake"="green3")) 
veg_median_greenness
ggsave("D:/SERDP/Mission_Creek/output_imagery/median_veg_greenness.png", veg_median_greenness, width=6, height=3, device="png")
# Deciduousness (from GEE Phenology)
# Max Greenness (from GEE Phenology)
veg_max_greenness <- ggplot(sample_topography_mc_rc) + 
  geom_smooth(aes(x=dist_up_creek,y=max_greenness,group=str_to_title(creek_name),col=str_to_title(creek_name)), level=0.95, method="loess") + 
  geom_point(aes(x=dist_up_creek,y=max_greenness,group=str_to_title(creek_name),col=str_to_title(creek_name)), alpha=0.5) + 
  geom_vline(data=geologic_boundaries, aes(xintercept=(position_rattlesnake+position_mission)/2), col="gray", linetype="dashed") + 
  geom_vline(data=mission_tunnel, aes(xintercept=position), col="magenta", alpha=0.5) +
  geom_vline(data=summer_water_boundaries, aes(xintercept=position, col=creek), alpha=0.5, linetype="dashed") +
  ggtitle("Maximum Vegetation Greenness by Longitudinal Position") + 
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  scale_x_continuous(limits=c(0,8000), expand=c(0,0)) +
  xlab("Longitudinal Position (m)") + 
  ylab("Maximum Greenness") +
  scale_color_manual(name="",values = c("Mission"="orange","Rattlesnake"="green3")) 
veg_max_greenness
ggsave("D:/SERDP/Mission_Creek/output_imagery/max_veg_greenness.png", veg_max_greenness, width=6, height=3, device="png")
# Deciduousness (from GEE Phenology)


# Generate t test comparisons between rattlesnake and mission by formation and variable type

testVariableInFormation <- function(variable_name, formation_name)
{
  test_result <- t.test((sample_topography_mc_rc %>% filter(as.character(creek_name)=="rattlesnake",as.character(formation)==formation_name))[,variable_name],
                        (sample_topography_mc_rc %>% filter(as.character(creek_name)=="mission",as.character(formation)==formation_name))[,variable_name])
  print("")
  print(paste("Difference in ", variable_name, " between creeks for formation ", formation_name))
  print(test_result)
  result_as_data_frame <- data.frame(variable = variable_name,
                                     formation = formation_name, 
                                     p_value = test_result$p.value,
                                     t_stat = test_result$parameter,
                                     estimate_rattlesnake = (test_result$estimate)[1],
                                     estimate_mission = (test_result$estimate)[2],
                                     difference = (test_result$estimate)[1] - (test_result$estimate)[2],
                                     conf_int_low = (test_result$conf.int)[1],
                                     conf_int_high = (test_result$conf.int)[2])
  return(result_as_data_frame)
}

all_tests <- rbind(testVariableInFormation("mean_return_height","Alluvium"),
                   testVariableInFormation("mean_return_height","Sespe"),
                   testVariableInFormation("mean_return_height","Coldwater"),
                   testVariableInFormation("mean_return_height","Cozydell"),
                   testVariableInFormation("mean_return_height","Matilija"),
                   testVariableInFormation("max_return_height","Alluvium"),
                   testVariableInFormation("max_return_height","Sespe"),
                   testVariableInFormation("max_return_height","Coldwater"),
                   testVariableInFormation("max_return_height","Cozydell"),
                   testVariableInFormation("max_return_height","Matilija"),
                   testVariableInFormation("decid_frac","Alluvium"),
                   testVariableInFormation("decid_frac","Sespe"),
                   testVariableInFormation("decid_frac","Coldwater"),
                   testVariableInFormation("decid_frac","Cozydell"),
                   testVariableInFormation("decid_frac","Matilija"),
                   testVariableInFormation("median_greenness","Alluvium"),
                   testVariableInFormation("median_greenness","Sespe"),
                   testVariableInFormation("median_greenness","Coldwater"),
                   testVariableInFormation("median_greenness","Cozydell"),
                   testVariableInFormation("median_greenness","Matilija"),
                   testVariableInFormation("max_greenness","Alluvium"),
                   testVariableInFormation("max_greenness","Sespe"),
                   testVariableInFormation("max_greenness","Coldwater"),
                   testVariableInFormation("max_greenness","Cozydell"),
                   testVariableInFormation("max_greenness","Matilija"))
num_wet_dry_by_formation <- sample_topography_mc_rc %>% 
  group_by(creek_name, formation, summer_dry) %>% 
  tally()


# code to look up checklist numbers and compare vegetation to bird records from audio




