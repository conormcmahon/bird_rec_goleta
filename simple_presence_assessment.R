
library(tidyverse)
library(lubridate)
library(sf)
library(unmarked)
library(janitor)
library(AICcmodavg)

# Detection Parameters
#   Confidence threshold - minimum confidence for a call to be counted
confidence_threshold <- 0.6
#   Detections Threshold - minimum number of detections for a bird to be flagged 'present' at a site
detections_threshold <- 100
# To be present, a bird must have at least one detection > confidence_threshold, and a total number of detections > detections_threshold

# Load detections data, including confidence value
all_detections <- rbind(read_csv("G:/Bioacoustics/Goleta_2022/detections/all_detections_compiled.csv"),
                        read_csv("G:/Bioacoustics/Goleta_2023/detections/all_detections_compiled_2023.csv"))

# Load survey effort data
survey_effort <- rbind(read_csv("G:/Bioacoustics/Goleta_2023/survey_time.csv"),
                       read_csv("G:/Bioacoustics/Goleta_2022/survey_time.csv")) %>%
  filter(year %in% c(2022, 2023))

# For one site where we experimentally tested several units, keep only one recorder type
#   First for detections
all_detections <- all_detections %>% 
  filter(!(site %in% c("SR_12_SM", "SR_12_AR")))
all_detections[all_detections$site == "SR_12_AM",]$site = "SR_12"
#   Next for survey effort
survey_effort <- survey_effort %>% 
  filter(!(site %in% c("SR_12_SM", "SR_12_AR")))
survey_effort[survey_effort$site == "SR_12_AM",]$site = "SR_12"


# Get hourly effort
hourly_effort <- survey_effort %>% 
  group_by(site, year, month, day, hour) %>%
  summarize(minutes = n()) %>%
  mutate(sitetime = paste(site, year, month, day, hour, sep="_")) 

# For each species, get the number of detections by hour
# Also, get number of detections with confidence > threshold
hourly_detections <- all_detections %>% 
  group_by(common_name, site, year, month, day, hour) %>% 
  summarize(num_detections = n(),
            num_good_detections = sum(confidence > confidence_threshold),
            max_conf = max(confidence),
            median_conf = median(confidence))

# For each hour and site surveyed, get number of detections of target species
getNormalizedCallRates <- function(target_species)
{
  # Filter to target species
  target_data <- hourly_detections %>% 
    filter(common_name == target_species)
  # Normalize number of detections / hr by number of minutes recorded
  detections_and_effort <- merge(hourly_effort, target_data, all=TRUE) %>% 
    mutate(num_detections_norm = num_detections / minutes,
           num_good_detections_norm = num_good_detections / minutes,
           doy = lubridate::yday(paste(year, sprintf("%02d", month), sprintf("%02d", day), sep="-")))
  detections_and_effort <- replace_na(detections_and_effort,
                                      list(common_name = target_species,
                                           num_detections = 0,
                                           num_good_detections = 0,
                                           max_conf = 0,
                                           median_conf = 0,
                                           num_detections_norm = 0,
                                           num_good_detections_norm = 0)) 
  return(detections_and_effort)
}


# Generate habitat data

# Surface water!
field_data <- st_read("D:/birdrec/field_data_goleta_2022/goleta_qfield_final_2023/birdrec_site_goleta_2022.gpkg") %>%
  drop_na(timestamp_collected) %>%                             # remove stolen or missing data 
  mutate(year = as.numeric(substr(timestamp_collected, 1, 4)), # add year
         site = str_replace(site_name, "-", "_")) %>%          # clean up formatting
  dplyr::select(-site_name)                                    # clean up formatting 
# Replace NA with FALSE
field_data$surface_flow <- replace_na(field_data$surface_flow, FALSE)
field_data$pooling_water <- replace_na(field_data$pooling_water, FALSE)
field_data$seeps_or_springs <- replace_na(field_data$seeps_or_springs, FALSE)
# Combine types of water into one surface water class
field_data$wet <- (field_data$surface_flow + field_data$pooling_water + field_data$seeps_or_springs) > 0

# Get variable to show whether water is ALWAYS present
field_data <- merge(as.data.frame(field_data), 
                    as.data.frame(field_data) %>% 
                      dplyr::select(-geometry) %>% 
                      group_by(site) %>% 
                      summarize(perennially_wet = sum(wet)==n()))

# Get distance to nearest water for dry sites 
getWaterDistance <- function(ind)
{
  # If target point itself was wet, set distance to 0
  if(field_data[ind,]$wet)
    return(0)
  # Filter out list of wet sites in this year
  wet_sites <- field_data %>% 
    filter(year == field_data[ind,]$year,
           wet) %>% 
    st_as_sf()
  # Get nearest wet point to target point
  nearest_wet_point_ind <- st_nearest_feature(field_data[ind,] %>% 
                                                st_as_sf(), 
                                              wet_sites)
  nearest_wet_point <- wet_sites[nearest_wet_point_ind,]
  # Get distance to nearest wet point
  return(as.numeric(st_distance(nearest_wet_point, 
                                field_data[ind,] %>% 
                                  st_as_sf())))
}
field_data$water_distance <- unlist(lapply(1:nrow(field_data), getWaterDistance))

surface_water <- as.data.frame(field_data) %>% 
  group_by(site, year) %>% 
  summarize(wet = sum(wet)>0,
            perennially_wet = sum(perennially_wet)==n(), 
            water_distance = max(water_distance))


# Remote Sensing
#   Loading this data from where it was constructed in ./remote_sensing_sampler.R
#   Come back and clean this up a lot later!
rs_data <- read_csv("D:/birdrec/remote_sensing_data/site_rs_data.csv") %>%
  mutate(site = site_name) %>%
  dplyr::select(-site_name)



# Point Counts
point_count_birds <- read_csv(here::here("summer_point_count_data.csv"))


# Get effort-normalized vocal frequency data on a target species
modelSpecies <- function(target_species, 
                         model_formula, 
                         starting_coefs, 
                         years_included=c(2022, 2023), 
                         hour_offset=2, 
                         use_point_counts=TRUE,
                         model_presence=TRUE,
                         quiet=FALSE)
{
  
  # Point Count Info
  point_count_observations <- point_count_birds %>% 
    mutate(year = as.numeric(substr(observation_date, 1, 4))) %>%
    filter(common_name == target_species,
           year %in% years_included) %>%
    mutate(site = str_replace(site, "-", "_"))
  # Get list of field site + year combinations where a bird was seen/heard in person on point counts
  field_validated_siteyears <- unique(paste(point_count_observations$site, 
                                            point_count_observations$year, sep="_"))
  
  # Generate call frequency data normalized by sampling time
  normalized_data <- getNormalizedCallRates(target_species) %>% 
    filter(year %in% years_included)
  
  # Get list of sites where the birds were detected by BirdNET
  sites_present <- normalized_data %>% 
    group_by(site, year) %>% 
    summarize(max_conf = max(max_conf), 
              num_detections = sum(num_detections),
              mean_detections = mean(num_detections_norm)) %>%
    filter(max_conf > confidence_threshold, num_detections > detections_threshold)
  # Optionally, add eBird sites (where birds were detected on point counts, but not by recorders)
  sites_detected_point_count_only <- point_count_observations %>%
    group_by(site, year) %>% 
    summarize(max_conf = 0,
              num_detections = 0,
              mean_detections = 0) %>%
    filter(!(paste(site,year) %in% paste(sites_present$site,sites_present$year)))
  # Optionally, update list of sites present with point count data
  if(use_point_counts)
    sites_present <- rbind(sites_present, sites_detected_point_count_only)
  if(!quiet)
    sites_present %>% View()
  
  normalized_data_with_covariates <- merge(merge(normalized_data, rs_data),
                                           as.data.frame(field_data) %>% 
                                             group_by(site, year) %>% 
                                             summarize(wet = sum(wet)>0,
                                                       surface_flow = sum(surface_flow)>0,
                                                       pooling_water = sum(pooling_water)>0,
                                                       perennially_wet = sum(perennially_wet)==n(),
                                                       water_distance = max(water_distance))) %>%
    filter(!(substr(site, 1, 2) %in% c("CS", "KR"))) %>%
    mutate(present = ((num_good_detections>0) *  # there were at least some detections
                        (paste(site,year) %in% paste(sites_present$site,sites_present$year))), # > 50 detections at this site total
           doy_norm = (doy-min(doy))/(max(doy)-min(doy)),
           year_norm = (year-min(year))/(max(year)-min(year)),
           creek=substr(site,1,2),
           hour_sinusoid = sin(((hour+hour_offset)%%24)/24*pi))
  
  # Split sites by year and treat as if they are independent observations
  normalized_data_with_covariates <- rbind(normalized_data_with_covariates %>% filter(year == 2022),
                                           normalized_data_with_covariates %>% filter(year == 2023))
  
  ggplot(normalized_data_with_covariates) + 
    geom_smooth(aes(x=hour, y=num_good_detections_norm)) + 
    facet_wrap(~(wet+0)+(doy>190))
  
  # Generate species presence matric
  presence <- normalized_data_with_covariates %>% 
    mutate(observation_period = paste("obs", month, day, hour, sep="_"),
           site_year = paste(site, year, sep="_")) %>% 
    #filter(year == 2022) %>%
    dplyr::select(site_year, observation_period, present) %>%
    pivot_wider(names_from=observation_period, values_from=present)
  if(!quiet)
    presence %>% View("Presence")
  
  # Habitat covariates
  site_covariates <- normalized_data_with_covariates %>% 
    mutate(observation_period = paste("obs", month, day, hour, sep="_"),
           site_year = paste(site, year, sep="_")) %>% 
    #filter(year == 2022) %>%
    #dplyr::select(c(16:57, ncol(normalized_data_with_covariates)+1, ncol(normalized_data_with_covariates)+2)) %>%
    group_by(site_year) %>%
    summarize(cluster_1 = mean(cluster_1, na.rm=TRUE),
              cluster_2 = mean(cluster_2, na.rm=TRUE),
              cluster_3 = mean(cluster_3, na.rm=TRUE),
              cluster_4 = mean(cluster_4, na.rm=TRUE),
              cluster_5 = mean(cluster_5, na.rm=TRUE),
              cluster_6 = mean(cluster_6, na.rm=TRUE),
              cluster_7 = mean(cluster_7, na.rm=TRUE),
              cluster_8 = mean(cluster_8, na.rm=TRUE),
              cluster_9 = mean(cluster_9, na.rm=TRUE),
              cluster_10 = mean(cluster_10, na.rm=TRUE),
              low_veg = mean(low_veg, na.rm=TRUE),
              mid_veg = mean(mid_veg, na.rm=TRUE),
              tall_veg = mean(tall_veg, na.rm=TRUE),
              max_height = mean(max_height, na.rm=TRUE),
              mean_height = mean(mean_height, na.rm=TRUE),
              mean_max_height = mean(mean_max_height, na.rm=TRUE),
              std_max_height = mean(std_max_height, na.rm=TRUE),
              decid = mean(decid, na.rm=TRUE),
              summer_greenness = mean(summer_greenness, na.rm=TRUE),
              dense_veg_fraction = mean(dense_veg_fraction, na.rm=TRUE),
              std_max_greenness = mean(std_max_greenness, na.rm=TRUE),
              mean_std_greenness = mean(mean_std_greenness, na.rm=TRUE),
              elevation = mean(elevation, na.rm=TRUE),
              std_elevation = mean(std_elevation, na.rm=TRUE),
              mean_slope = mean(mean_slope),
              max_slope = mean(max_slope),
              min_slope = mean(min_slope),
              std_slope = mean(std_slope),
              height_pct_95 = mean(height_pct_95, na.rm=TRUE),
              height_pct_90 = mean(height_pct_90, na.rm=TRUE),
              height_pct_80 = mean(height_pct_80, na.rm=TRUE),
              height_pct_70 = mean(height_pct_70, na.rm=TRUE),
              height_pct_60 = mean(height_pct_60, na.rm=TRUE),
              height_pct_50 = mean(height_pct_50, na.rm=TRUE),
              height_pct_40 = mean(height_pct_40, na.rm=TRUE),
              height_pct_30 = mean(height_pct_30, na.rm=TRUE),
              height_pct_20 = mean(height_pct_20, na.rm=TRUE),
              height_pct_10 = mean(height_pct_10, na.rm=TRUE),
              wet = sum(wet) > 0,
              perennially_wet = sum(perennially_wet)==n(), 
              surface_flow = sum(surface_flow) > 0,
              pooling_water = sum(pooling_water) > 0,
              water_distance = max(water_distance),
              site = modal(site),
              creek = modal(creek)) 
  # Normalize predictor variables
  site_covariates <- site_covariates %>%
    mutate(max_height = max_height/max(max_height),
           mean_height = mean_height/max(mean_height),
           mean_max_height = mean_max_height/max(mean_max_height),
           std_max_height = std_max_height/max(std_max_height),
           elevation = elevation/max(elevation),
           std_elevation = std_elevation/max(std_elevation),
           mean_slope = mean_slope/max(mean_slope),
           max_slope = max_slope/max(max_slope),
           min_slope = min_slope/max(min_slope),
           std_slope = std_slope/max(std_slope),
           height_pct_ratio_50_90 = height_pct_50 / height_pct_90,
           water_distance = water_distance / max(water_distance))

  # Observation data
  #  Hour of each observation block
  observation_hour <- normalized_data_with_covariates %>% 
    mutate(observation_period = paste("obs", month, day, hour, sep="_"),
           site_year = paste(site, year, sep="_")) %>% 
    #filter(year == 2022) %>%
    dplyr::select(site_year, observation_period, hour) %>%
    pivot_wider(names_from=observation_period, values_from=hour) %>%
    dplyr::select(-site_year)
  #     Normalize observation hour to fit a sinusoid with peak at 10am
  observation_hour_normalized <- sin(((observation_hour+hour_offset)%%24)/24*pi)
  #  Day of Year for each observation block (0 to 365)
  observation_doy <- normalized_data_with_covariates %>% 
    mutate(observation_period = paste("obs", month, day, hour, sep="_"),
           site_year = paste(site, year, sep="_")) %>% 
    #filter(year == 2022) %>%
    dplyr::select(site_year, observation_period, doy) %>%
    pivot_wider(names_from=observation_period, values_from=doy) %>%
    dplyr::select(-site_year)
  #     Normalize observation DOY to range from 0 to 1
  observation_doy_normalized <- (observation_doy-min(observation_doy,na.rm=TRUE)) / (max(observation_doy, na.rm=TRUE)-min(observation_doy,na.rm=TRUE))
  
  # Build occupancy input dataset
  occupancy_data<- unmarkedFrameOccu(y=as.matrix(presence[,2:ncol(presence)]),
                                     obsCovs = list(hour = observation_hour_normalized,
                                                    doy = observation_doy_normalized),
                                     siteCovs = site_covariates)
  
  if(use_point_counts)
  {
    # Get index of sites where birds were detected on point counts, within list of rows in site_covariates
    field_validated_site_indices <- which(site_covariates$site_year %in% field_validated_siteyears)
    if(!quiet)
      print(paste("  The following point counts had this species: ", 
                paste(site_covariates[field_validated_site_indices,]$site_year,collapse=" "), sep=""))
    
    # Run occupancy model 
    if(model_presence)
      occupancy_model <- occu(formula = model_formula, 
                              data = occupancy_data,
                              starts = starting_coefs,
                              knownOcc = field_validated_site_indices)
    else
      occupancy_model <- occu(formula = model_formula, 
                              data = occupancy_data,
                              starts = starting_coefs,
                              knownOcc = field_validated_site_indices,
                              linkPsi = 'cloglog')
  }
  else
  {
    # Run occupancy model 
    if(model_presence)
      occupancy_model <- occu(formula = model_formula, 
                              data = occupancy_data,
                              starts = starting_coefs)
    else
      occupancy_model <- occu(formula = model_formula, 
                              data = occupancy_data,
                              starts = starting_coefs,
                              linkPsi = 'cloglog')
  }
  if(!quiet)
    print(summary(occupancy_model))
  sites_detected <- sum(rowSums(presence[,2:ncol(presence)], na.rm=TRUE) > 0)
  sites_detected_frac <- sum(rowSums(presence[,2:ncol(presence)], na.rm=TRUE) > 0)/nrow(presence)
  if(!quiet)
    print(paste("This bird was detected at ", sites_detected, " out of ", nrow(presence), " sites, or ", round(sites_detected_frac*100,1), "%", sep=""))
  
  
  return(list(occ_mod = occupancy_model,
              occ_data = occupancy_data,
              norm_data_with_covar = normalized_data_with_covariates,
              site_covariates = site_covariates,
              present_list = sites_present,
              presence_mat = presence[,2:ncol(presence)],
              known_sites = field_validated_site_indices))
}


testModelsForSpecies <- function(target_species)
{
  model_all <- modelSpecies(target_species,
                            ~hour + doy  # detection model
                            ~wet + decid + summer_greenness + height_pct_90 + std_max_height + low_veg,
                            rep(-1, 10))
}

# TO DO IN MORNING
# Finish reading tutorial --> look for caveats 
# XXXXX     Data Prep
# XXXXX    Gather LiDAR data
# XXXXX      Gather phenology data
# XXXXX      Gather veg field survey data
# XXXXX      Gather elevation, creek gradient?? 
# XXXXX     Streamline model selection?
# XXXXX       Does it always make sense to do... 
# XXXXX         Presence:
# XXXXX           willow_frac, oak_frac, eucalyptus_frac?   --> no, we're missing this at lots of sites
# XXXXX          understory density, canopy height, canopy roughness?
# XXXXX          surface water presence
# XXXXX           (berry plants)
# XXXXX           (phenology patterns)
# XXXXX           (LST or PRISM --> seasonal climatology)
# XXXXX         Detection:
# XXXXX           hour (sinusoid)
# XXXXX           day of year
# XXXXX           flowing water * creek gradient? 
# XXXXX           (temperature / cloud cover --> met station)
# XXXXX     Streamline best-hour selection? 
# XXXXX     Iterate on confidence thresholds? 
# XXXXX     Lock down a set of models for this poster! 
# XXXXX     Start generating some figures 
# XXXXX     Map of predicted habitat????


newmod4 <- modelSpecies("Black-headed Grosbeak",
                        ~hour + doy + doy:wet  # detection model
                        ~wet + decid,
                        rep(-1, 7),
                        years_included = 2023,
                        use_point_counts = TRUE)

test_df <- data.frame(hour = rep(sin((((1:50)/50+(2/24))%%1)*pi),2),
                      doy = c(rep(0,50),rep(1,50)),
                      wet = rep(0,100),
                      hour_real = rep(1:50/50*24,2))
predictions <- modavgPred(list(newmod4[[1]]),
                          parm.type = "detect",
                          newdata = test_df)
prediction_df <- as.data.frame(predictions)
prediction_df <- cbind(prediction_df, test_df)
max(prediction_df$mod.avg.pred)
mean(prediction_df$mod.avg.pred)
median(prediction_df$mod.avg.pred)
ggplot(prediction_df) +
  geom_line(aes(x=hour_real, y=lower.CL, group=doy, linetype=(doy==1)), col="red") +
  geom_line(aes(x=hour_real, y=mod.avg.pred, group=doy, linetype=(doy==1)), col="black") +
  geom_line(aes(x=hour_real, y=upper.CL, group=doy, linetype=(doy==1)), col="blue")

test_df <- data.frame(wet = c(rep(0,50),rep(1,50)),
decid = rep((1:50)/50,2),
low_veg = rep(.5,100))



