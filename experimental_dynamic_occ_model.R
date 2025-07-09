


# Get effort-normalized vocal frequency data on a target species
modelSpeciesDynamic <- function(target_species, 
                                occupancy_formula,
                                colonization_formula,
                                extinction_formula,
                                detection_formula, 
                                starting_coefs, 
                                years_included=c(2022, 2023), 
                                hour_offset=2, 
                                use_point_counts=TRUE,
                                model_presence=TRUE)
{
  
  # Point Count Info
  point_count_observations <- point_count_birds %>% 
    mutate(year = as.numeric(substr(observation_date, 1, 4))) %>%
    filter(common_name == target_species,
           year == min(years_included)) %>%
    mutate(site = str_replace(site, "-", "_"))
  
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
  sites_present %>% View("Sites Present")
  
  normalized_data_with_covariates <- merge(merge(normalized_data, rs_data, all=TRUE),
                                           as.data.frame(field_data) %>% 
                                             group_by(site, year) %>% 
                                             summarize(wet = sum(wet)>0,
                                                       surface_flow = sum(surface_flow)>0,
                                                       pooling_water = sum(pooling_water)>0,
                                                       perennially_wet = sum(perennially_wet)==n(),
                                                       water_distance = max(water_distance)), 
                                           all=TRUE) %>%
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
    mutate(observation_period = paste("obs", year, month, day, hour, sep="_")) %>% 
    dplyr::select(site, observation_period, present) %>%
    pivot_wider(names_from=observation_period, values_from=present)
  presence %>% View("Presence")
  
  # Habitat covariates
  site_covariates <- normalized_data_with_covariates %>% 
    mutate(observation_period = paste("obs", year, month, day, hour, sep="_")) %>% 
    #dplyr::select(c(16:57, ncol(normalized_data_with_covariates)+1, ncol(normalized_data_with_covariates)+2)) %>%
    group_by(site) %>%
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
              perennially_wet = sum(perennially_wet)==n(),
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
           height_pct_ratio_50_90 = height_pct_50 / height_pct_90)
  
  # Site-Year Covariates - Water
  #   Year
  yearly_year <- data.frame(site = unique(normalized_data_with_covariates$site))
  yearly_year$x1 <- 2022
  yearly_year$x2 <- 2023
  
  #   Surface water (any)
  yearly_wet <- normalized_data_with_covariates %>%
    group_by(site, year) %>%
    summarize(wet = max(wet) > 0) %>%
    pivot_wider(names_from=year, values_from=wet)
  #   Surface Flow (active flow in channel)
  yearly_flow <- normalized_data_with_covariates %>%
    group_by(site, year) %>%
    summarize(perennially_wet = max(perennially_wet)==n()) %>%
    pivot_wider(names_from=year, values_from=perennially_wet)
  #   Pooling Water (pools but no active flow)
  yearly_pools <- normalized_data_with_covariates %>%
    group_by(site, year) %>%
    summarize(pooling_water = sum(pooling_water) > 0) %>%
    pivot_wider(names_from=year, values_from=pooling_water)
  #   Distance to Surface Water (for dry sites)
  yearly_water_dist <- normalized_data_with_covariates %>%
    group_by(site, year) %>%
    summarize(water_distance = max(water_distance)) %>%
    pivot_wider(names_from=year, values_from=water_distance)
  #   Normalize water distances by maximum in whole dataset
  yearly_water_dist[,2:ncol(yearly_water_dist)] <- yearly_water_dist[,2:ncol(yearly_water_dist)] / 
    max(yearly_water_dist[,2:ncol(yearly_water_dist)], na.rm=TRUE)
  
  yearly_wet[,1:ncol(yearly_wet)] %>% View()
  
  # Observation data
  #  Hour of each observation block
  observation_hour <- normalized_data_with_covariates %>% 
    mutate(observation_period = paste("obs", year, month, day, hour, sep="_")) %>% 
    dplyr::select(site, observation_period, hour) %>%
    pivot_wider(names_from=observation_period, values_from=hour) %>%
    dplyr::select(-site)
  #     Normalize observation hour to fit a sinusoid with peak at 10am
  observation_hour_normalized <- sin(((observation_hour+hour_offset)%%24)/24*pi)
  #  Day of Year for each observation block (0 to 365)
  observation_doy <- normalized_data_with_covariates %>% 
    mutate(observation_period = paste("obs", year, month, day, hour, sep="_")) %>% 
    dplyr::select(site, observation_period, doy) %>%
    pivot_wider(names_from=observation_period, values_from=doy) %>%
    dplyr::select(-site)
  #     Normalize observation DOY to range from 0 to 1
  observation_doy_normalized <- (observation_doy-min(observation_doy,na.rm=TRUE)) / (max(observation_doy, na.rm=TRUE)-min(observation_doy,na.rm=TRUE))
  
  
  # Correct number of observations to be equivalent both years (pad with NA)
  #  First, in the presence dataset
  presence_2022 <- presence[,grep("2022", names(presence))]
  presence_2023 <- presence[,grep("2023", names(presence))]
  presence <- cbind(presence[,1], 
                    presence_2022,
                    matrix(rep(NA,(ncol(presence_2023)-ncol(presence_2022))*nrow(presence)), nrow=nrow(presence)),
                    presence_2023)
  names(presence) <- c("site", 
                       paste("obs_2022_", 1:ncol(presence_2023), sep=""),
                       paste("obs_2023_", 1:ncol(presence_2023), sep=""))
  #  Next, the same for observation hour
  hour_2022 <- observation_hour_normalized[,grep("2022", names(observation_hour_normalized))]
  hour_2023 <- observation_hour_normalized[,grep("2023", names(observation_hour_normalized))]
  observation_hour_normalized <- cbind(hour_2022,
                                       matrix(rep(NA,(ncol(hour_2023)-ncol(hour_2022))*nrow(observation_hour_normalized)), nrow=nrow(observation_hour_normalized)),
                                       hour_2023)
  names(observation_hour_normalized) <- c(paste("obs_2022_", 1:ncol(hour_2023), sep=""),
                                          paste("obs_2023_", 1:ncol(hour_2023), sep=""))
  #  Next, the same for observation day-of-year
  doy_2022 <- observation_doy_normalized[,grep("2022", names(observation_doy_normalized))]
  doy_2023 <- observation_doy_normalized[,grep("2023", names(observation_doy_normalized))]
  observation_doy_normalized <- cbind(doy_2022,
                                      matrix(rep(NA,(ncol(doy_2023)-ncol(doy_2022))*nrow(observation_doy_normalized)), nrow=nrow(observation_doy_normalized)),
                                      doy_2023)
  names(observation_doy_normalized) <- c(paste("obs_2022_", 1:ncol(doy_2023), sep=""),
                                         paste("obs_2023_", 1:ncol(doy_2023), sep=""))
  
  # Build occupancy input dataset
  occupancy_data<- unmarkedMultFrame(y=as.matrix(presence[,2:ncol(presence)]),
                                     obsCovs = list(hour = observation_hour_normalized,
                                                    doy = observation_doy_normalized),
                                     siteCovs = site_covariates,
                                     yearlySiteCovs = data.frame(site = (yearly_year %>% pivot_longer(2:3, names_to="var", values_to="year") %>% arrange(site,var) %>% dplyr::select(-var))$site,
                                                                 year = (yearly_year %>% pivot_longer(2:3, names_to="var", values_to="year") %>% arrange(site,var) %>% dplyr::select(-var))$year,
                                                                 wet = (yearly_wet %>% pivot_longer(2:3, names_to="var", values_to="wet") %>% arrange(site,var) %>% dplyr::select(-var))$wet,
                                                                 surface_flow = (yearly_flow %>% pivot_longer(2:3, names_to="var", values_to="surface_flow") %>% arrange(site,var) %>% dplyr::select(-var))$surface_flow,
                                                                 pooling = (yearly_pools %>% pivot_longer(2:3, names_to="var", values_to="pooling") %>% arrange(site,var) %>% dplyr::select(-var))$pooling,
                                                                 distance_to_water = (yearly_water_dist %>% pivot_longer(2:3, names_to="var", values_to="water_dist") %>% arrange(site,var) %>% dplyr::select(-var))$water_dist),
                                     numPrimary = 2)
# 
#   return(list(presence,
#               observation_hour_normalized,
#               observation_doy_normalized,
#               site_covariates,
#               data.frame(site = (yearly_year %>% pivot_longer(2:3, names_to="var", values_to="year") %>% arrange(site,var) %>% dplyr::select(-var))$site,
#                          year = (yearly_year %>% pivot_longer(2:3, names_to="var", values_to="year") %>% arrange(site,var) %>% dplyr::select(-var))$year,
#                          wet = (yearly_wet %>% pivot_longer(2:3, names_to="var", values_to="wet") %>% arrange(site,var) %>% dplyr::select(-var))$wet,
#                          surface_flow = (yearly_flow %>% pivot_longer(2:3, names_to="var", values_to="surface_flow") %>% arrange(site,var) %>% dplyr::select(-var))$surface_flow,
#                          pooling = (yearly_pools %>% pivot_longer(2:3, names_to="var", values_to="pooling") %>% arrange(site,var) %>% dplyr::select(-var))$pooling,
#                          distance_to_water = (yearly_water_dist %>% pivot_longer(2:3, names_to="var", values_to="water_dist") %>% arrange(site,var) %>% dplyr::select(-var))$water_dist),
#               normalized_data_with_covariates,
#               occupancy_data))
#   
  if(use_point_counts)
  {
    # Get index of sites where birds were detected on point counts, within list of rows in site_covariates
    field_validated_site_indices <- which(site_covariates$site %in% point_count_observations$site)
    print(paste("  The following point counts had this species in year 1: ", 
                paste(site_covariates[field_validated_site_indices,]$site,collapse=" "), sep=""))
    
    # Run occupancy model 
    if(model_presence)
    {
      occupancy_model <- colext(psiformula = occupancy_formula,
                                gammaformula = colonization_formula,
                                epsilonformula = extinction_formula,
                                pformula = detection_formula,
                                data = occupancy_data,
                                starts = starting_coefs,
                                knownOcc = field_validated_site_indices)
    }
    else
    {
      occupancy_model <- colext(psiformula = occupancy_formula,
                                gammaformula = colonization_formula,
                                epsilonformula = extinction_formula,
                                pformula = detection_formula,
                                data = occupancy_data,
                                starts = starting_coefs,
                                knownOcc = field_validated_site_indices,
                                linkPsi = 'cloglog')
    }
  }
  else
  {
    # Run occupancy model 
    if(model_presence)
      occupancy_model <- colext(psiformula = occupancy_formula,
                                gammaformula = colonization_formula,
                                epsilonformula = extinction_formula,
                                pformula = detection_formula,
                                data = occupancy_data,
                                starts = starting_coefs)
    else
      occupancy_model <- colext(psiformula = occupancy_formula,
                                gammaformula = colonization_formula,
                                epsilonformula = extinction_formula,
                                pformula = detection_formula,
                                data = occupancy_data,
                                starts = starting_coefs,
                                linkPsi = 'cloglog')
  }
  print(summary(occupancy_model))
  sites_detected <- sum(rowSums(presence[,2:ncol(presence)], na.rm=TRUE) > 0)
  sites_detected_frac <- sum(rowSums(presence[,2:ncol(presence)], na.rm=TRUE) > 0)/nrow(presence)
  print(paste("This bird was detected at ", sites_detected, " out of ", nrow(presence), " sites, or ", round(sites_detected_frac*100,1), "%", sep=""))
  
  
  return(list(occupancy_model,
              occupancy_data,
              normalized_data_with_covariates,
              site_covariates,
              sites_present,
              presence[,2:ncol(presence)]))
}


temp <- modelSpeciesDynamic("Black-headed Grosbeak",
                            ~1,
                            ~1,
                            ~1,
                            ~1,
                            rep(1, 4),
                            use_point_counts = TRUE)

