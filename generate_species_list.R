

# Generate species list for a given site

# Specify these parameters!
#target_site <- c("MC")
#output_name <- "botanic_garden"
target_site <- c("AC","BR","EW","LL","MC","NC","RC","SA","SJ","SR")
output_name <- "all_data"
confidence_threshold <- 0.85
detections_threshold <- 250

# Required libraries
library(tidyverse)
library(sf)

# Load list of bird names
species_metadata <- read_csv("D:/birdrec/scripts/birdnet_species_list_processing/all_species_metadata.csv")
bird_names <- (species_metadata %>% filter(bird))$common_name

# Load a list which includes the species from the Riparian Broad species list
birdnet_riparian_broad_list <- read_csv("D:/birdrec/protocol/birdnet_riparian_broad.csv")

# Load another list which includes commonly-reported birds from eBird in the county which are missing above
missing_from_birdnet_list <- read_csv("D:/birdrec/protocol/missing_from_birdnet_riparian_broad.csv") %>%
  filter(num_reports > 2500)

# Combine the above into an 'allowed species' list
allowed_species <- unique(c(birdnet_riparian_broad_list$common_name,
                            missing_from_birdnet_list$common_name,
                            "Mountain Quail",
                            "Greater Roadrunner",
                            "Cassin's Vireo",
                            "Olive-sided Flycatcher"))

# Load data on survey effort (number of hours recorded)
effort_22 <- read_csv("G:/Bioacoustics/Goleta_2022/survey_time.csv")
effort_23 <- read_csv("G:/Bioacoustics/Goleta_2023/survey_time.csv")
effort_24 <- read_csv("G:/Bioacoustics/Goleta_2024/survey_time.csv")
effort_overall <- rbind(effort_22, effort_23, effort_24)
effort_summary <- effort_overall %>%
  group_by(site) %>%
  summarize(survey_days = n()/60/24)

# Load spatial data on sites
field_sites <- st_read("D:/birdrec/field_data_goleta_2022/goleta_qfield_final_2023/birdrec_site_goleta_2022.gpkg") %>%
  mutate(site = site_name) %>%
  select(-site_name)
field_sites$site <- str_replace(field_sites$site, "-", "_") 
field_sites_averaged <- field_sites %>% aggregate(., 
                                                  by = list(.$site), 
                                                  function(x) x = x[1]) %>% 
  st_centroid() %>% 
  select(-Group.1)

# Load point counts
point_counts <- read_csv("D:/birdrec/scripts/bird_rec_goleta/summer_point_count_data.csv") %>%
  filter(substr(site, 1,2) %in% target_site)
point_counts$distant <- grepl("distant", point_counts$species_comments)
point_count_species_list <- unique((point_counts %>% filter(substr(site, 1,2) %in% target_site))$common_name)

# Load detections

all_detections <- rbind(read_csv("G:/Bioacoustics/Goleta_2022/detections_compiled/all_detections.csv"),
                        read_csv("G:/Bioacoustics/Goleta_2023/detections_compiled/all_detections.csv"),
                        read_csv("G:/Bioacoustics/Goleta_2024/detections_compiled/all_detections.csv")) %>%
  filter(substr(site,1,2) %in% target_site,
         common_name %in% bird_names) %>%
  mutate(timestamp = paste(sprintf("%04.0f",year),
                           sprintf("%02.0f",month),
                           sprintf("%02.0f",day),
                           "_",
                           sprintf("%02.0f",hour),
                           sprintf("%02.0f",minute),
                           sprintf("%02.0f",second),
                           sep=""),
         site_broad = substr(site, 1,2)) %>%
  filter(common_name %in% allowed_species)
all_detections <- merge(all_detections, 
                        effort_summary, 
                        by="site")

# Get overall summary of ARU data across all sites
num_recording_sites <- length(unique(all_detections$site))
aru_summary <- all_detections %>%
  group_by(common_name) %>% 
  summarize(num_detections = n(),
            num_detections_high_confidence = sum(confidence > confidence_threshold),
            num_sites = length(unique(site)),
            fraction_of_sites = round(length(unique(site))/num_recording_sites,3),
            conf_q_50 = median(confidence),
            conf_q_90 = quantile(confidence, 0.9),
            conf_q_100 = max(confidence))
aru_summary_by_site <- all_detections %>%
  filter(confidence > confidence_threshold) %>%
  group_by(common_name, site) %>%
  summarize(detected = 1) %>%
  group_by(common_name) %>%
  summarize(num_sites_high_confidence = n(),
            fraction_of_sites_high_confidence = n()/num_recording_sites)
aru_summary <- merge(aru_summary, aru_summary_by_site,
                     by="common_name", all=TRUE)
# Check over initial list and look for suspect species - especially those with few detections or low confidence
View(aru_summary %>% 
       arrange(-num_detections))
# Look at lists of detections for individual suspect species, and check the records manually
all_detections %>% filter(common_name == "Bullock's Oriole") %>% 
  arrange(-confidence)

# From manual inspection for Ellwood Mesa, it looks like good cutoffs are 0.85 min(max(confidence)) and detections > 250
aru_species_list <- unique((aru_summary %>% filter(conf_q_100 > confidence_threshold, num_detections > detections_threshold))$common_name)
aru_summary <- aru_summary %>%
  filter(common_name %in% aru_species_list)

# Point Count summary
num_point_counts <- length(unique(point_counts$checklist_id))
point_count_summary <- point_counts %>%
  group_by(common_name) %>% 
  summarize(point_count_presences = n(),
            fraction_of_counts = round(n()/num_point_counts,3),
            total_individuals = sum(observation_count))
View(point_count_summary)

# Combine point count and ARU summaries by species (across all sites)
all_summaries <- merge(aru_summary, point_count_summary, by="common_name", all=TRUE) %>%
  select(c("common_name", 
           "num_detections_high_confidence", "num_sites_high_confidence", "fraction_of_sites_high_confidence",
           "conf_q_100", "conf_q_90", "conf_q_50",
           "total_individuals", "point_count_presences", "fraction_of_counts",
           "num_detections", "num_sites", "fraction_of_sites")) %>%
  arrange(-num_sites_high_confidence, 
          -num_detections_high_confidence, 
          -point_count_presences) %>%
  replace(is.na(.), 0)
View(all_summaries %>% arrange(-num_detections))
#write_csv(all_summaries, paste("D:/birdrec/reports/species_lists/creek_level/", output_name, ".csv", sep=""))

# Which species were picked up by recorders but not point counts?
aru_only_species <- aru_species_list[!(aru_species_list %in% point_count_species_list)]
print("Species only detected by ARUs:")
aru_only_species
# Which species were picked up by point counts but not recorders?
point_count_only_species <- point_count_species_list[!(point_count_species_list %in% aru_species_list)]
print("Species only detected by Point Counts:")
point_count_only_species

# Generate spatial summary information for a given species
getSpeciesSummary <- function(target_species, sites, count_threshold=0, confidence_threshold=0) 
{
  # Generate summary for each site
  target_data <- all_detections %>% 
    filter(substr(site,1,2) %in% sites, common_name == target_species) %>% 
    group_by(site) %>%
    summarize(count = n(),
              count_norm = n()/mean(survey_days),
              conf_q_100 = quantile(confidence, 1.0),
              conf_q_90 = quantile(confidence, 0.9),
              conf_q_80 = quantile(confidence, 0.8),
              conf_q_70 = quantile(confidence, 0.7),
              conf_q_60 = quantile(confidence, 0.6),
              conf_q_50 = quantile(confidence, 0.5),
              conf_q_40 = quantile(confidence, 0.4),
              conf_q_30 = quantile(confidence, 0.3),
              conf_q_20 = quantile(confidence, 0.2),
              conf_q_10 = quantile(confidence, 0.1)) %>%
    filter(count > count_threshold,
           conf_q_100 > confidence_threshold)
  # Merge with spatial data 
  target_data <- merge(target_data, field_sites_averaged)
  # Output to .gpkg file
  #st_write(target_data, paste("D:/birdrec/bird_maps/", target_species, ".gpkg", sep=""), append=FALSE)
  return(target_data)
}











