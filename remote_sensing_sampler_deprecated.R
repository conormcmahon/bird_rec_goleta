
# Code to sample remote sensing products for vegetation: 
#   Phenology (Sentinel-2 processed in GEE to 12-month greenness timeseries)
#   Structure (LiDAR height histograms)
# NOTE - need to come back in and add:
#   phenology data from 2023
#   LiDAR data from Cold Springs Creek (added in 2023)

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

# Load points shapefile
field_points <- st_read("D:/birdrec/field_data_goleta_2022/goleta_qfield_final_2023/birdrec_site_goleta_2022.gpkg") %>% 
  mutate(creek_name = substring(site_name,1,(nchar(as.character(site_name))-3))) %>%
  #drop_na(timestamp_deployed) %>% # remove sites which were point-count only (Kinevan Rd.) 
  mutate(site_name = str_replace(site_name, "-", "_"))
# Derive name for creek from site ID and add it to data frame
creek_strings <- unlist(lapply(str_split(field_points$site_name, "-"), function(str_list){return(str_list[[1]])}))
field_points$creek <- creek_strings

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
field_points_sym <-  field_points %>% filter(!(creek_name %in% c("BR","EW","NCOS","SR","AC","MYC"))) %>%
  filter(creek_name != "CSC")

# Sample rasters in radius around points
#   LiDAR
sample_radius <- 30
aq_samples <- raster::extract(readAll(stack(lidar_aq, crop(datastack_aq,lidar_aq))), field_points_aq, method='simple', buffer=sample_radius, cellnumbers=TRUE) #, fun=mean, na.rm=TRUE)
copr_samples <- raster::extract(readAll(stack(lidar_copr, crop(datastack_copr,lidar_copr))), field_points_copr, method='simple', buffer=sample_radius, cellnumbers=TRUE)
sr_samples <- raster::extract(readAll(stack(lidar_sr, crop(datastack_sr,lidar_sr))), field_points_sr, method='simple', buffer=sample_radius, cellnumbers=TRUE)
ac_myc_samples <- raster::extract(readAll(stack(lidar_ac_myc, crop(datastack_sym,lidar_ac_myc))), field_points_ac_myc, method='simple', buffer=sample_radius, cellnumbers=TRUE)
sym_samples <- raster::extract(readAll(stack(lidar_sym, crop(datastack_sym,lidar_sym))), field_points_sym, method='simple', buffer=sample_radius, cellnumbers=TRUE)

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
vegetation_df <- bind_rows(lapply(1:length(vegetation_samples), addSiteNameToSamples))

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

vegetation_structure_long <- vegetation_df %>% 
  dplyr::select(c(paste("bin_",1:30,sep=""),"site_name","creek_name","wet")) %>%
  pivot_longer(1:30,names_to="bin_string",values_to="fraction") %>%
  mutate(bin = as.numeric(substring(bin_string, 5,10)))
vegetation_structure_long <- merge(vegetation_structure_long,(sample_topography %>% dplyr::select(-creek_name)),by="site_name")
vegetation_phenology_long <- vegetation_df %>% 
  dplyr::select(c(paste("month_",1:12,sep=""),"site_name","creek_name","wet")) %>%
  pivot_longer(1:12,names_to="month_string",values_to="NDVI") %>%
  mutate(month = as.numeric(substring(month_string, 7,10)))
vegetation_phenology_long <- merge(vegetation_phenology_long,(sample_topography %>% dplyr::select(-creek_name)),by="site_name")

ggplot(vegetation_structure_long %>%
         filter(creek_name %in% c("RC","MC"))) + 
  geom_density_2d_filled(aes(y=bin, x=fraction), contour_var = "ndensity") + 
  facet_wrap(~site_name) + 
  scale_x_continuous(limits=c(.01,0.3)) + 
  scale_fill_manual(values=c(colorRampPalette(paletteer::paletteer_dynamic("cartography::green.pal", 
                                                                           n = 20, 
                                                                           direction = 1))(10)), 
                    aesthetics = c("fill", "color"))

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
            std_max_height = sd(max_height),
            mean_max_height = mean(max_height),
            max_height = max(max_height),
            mean_height = mean(mean_height),
            max_over_mean_height = max_height / mean_height,
            decid = sum(summer_ndvi > winter_ndvi)/n(),
            summer_greenness = sum(summer_ndvi*(summer_ndvi>0.4), na.rm=TRUE)/sum((summer_ndvi>0.4), na.rm=TRUE),
            dense_veg_fraction = sum(summer_ndvi > 0.4)/n(),
            std_max_greenness = sd(max_ndvi),
            mean_std_greenness = mean(stdev_ndvi),
            elevation = mean(dem),
            std_elevation = sd(dem),
            mean_slope = mean(slope),
            max_slope = max(slope),
            min_slope = min(slope),
            std_slope = sd(slope)) 
veg_df_collapsed$height_pct_95 <- getHeightQuantile(0.95) # /veg_df_collapsed$max_height    Was originally normalizing by max height, 
veg_df_collapsed$height_pct_90 <- getHeightQuantile(0.9) # /veg_df_collapsed$max_height     but actually should probably do this on the basis
veg_df_collapsed$height_pct_80 <- getHeightQuantile(0.8) # /veg_df_collapsed$max_height     of OVERALL max height, not per-site
veg_df_collapsed$height_pct_70 <- getHeightQuantile(0.7) # /veg_df_collapsed$max_height
veg_df_collapsed$height_pct_60 <- getHeightQuantile(0.6) # /veg_df_collapsed$max_height
veg_df_collapsed$height_pct_50 <- getHeightQuantile(0.5) # /veg_df_collapsed$max_height
veg_df_collapsed$height_pct_40 <- getHeightQuantile(0.4) # /veg_df_collapsed$max_height
veg_df_collapsed$height_pct_30 <- getHeightQuantile(0.3) # /veg_df_collapsed$max_height
veg_df_collapsed$height_pct_20 <- getHeightQuantile(0.2) # /veg_df_collapsed$max_height
veg_df_collapsed$height_pct_10 <- getHeightQuantile(0.1) # /veg_df_collapsed$max_height
veg_df_collapsed$summer_greenness <- replace(veg_df_collapsed$summer_greenness, is.na(veg_df_collapsed$summer_greenness), 0)

write_csv(veg_df_collapsed, "D:/birdrec/remote_sensing_data/site_rs_data.csv")




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




