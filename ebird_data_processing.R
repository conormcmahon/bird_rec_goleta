
library(auk)
library(tidyverse)
library(here)
library(sf)

# path to the ebird data file, here all Santa Barbara County data up until July 2022
# get the path to the example data included in the package
f_in <- here::here("..","..","eBird","ebird_ca_dec_2023","ebd_US-CA_smp_relOct-2023.txt")
# output text file
f_out <- "all_sb_co_data.txt"
ebird_data <- f_in %>% 
  # 1. reference file
  auk_ebd() %>% 
  # 2. define filters
#  auk_species(species = "Oak Titmouse") %>% 
  auk_county("US-CA-083") %>%
  auk_date(date=c("2022-06-01","2023-08-15")) %>% 
  # 3. run filtering
  auk_filter(file = f_out, overwrite=TRUE) %>% 
  # 4. read text file into r data frame
  read_ebd()

# Load field data on sites
#   Correct point count checklist field name to match eBird data
field_sites <- as.data.frame(st_read("D:/birdrec/field_data_goleta_2022/goleta_qfield_final_2023/birdrec_site_goleta_2022.gpkg")) %>%
  mutate(checklist_id = point_count_checklist_id,
         site = site_name) %>%
  select(-point_count_checklist_id, -site_name) %>%
  drop_na(checklist_id)

# Get list of all checklists in field dataset
valid_checklists <- unique(field_sites$checklist_id)
                              

# Get eBird data from field checklists
summer_survey_data <- ebird_data %>% 
  filter(checklist_id %in% valid_checklists) 

# Merge eBird data with field information
summer_survey_data <- merge(summer_survey_data, 
                            field_sites %>% 
                              dplyr::select(checklist_id, site, point_count_timestamp, timestamp_deployed, timestamp_collected))

# Write eBird and field data combination  to disk
write_csv(summer_survey_data, "summer_point_count_data.csv")


# 
# ebird_and_field_data <- merge(summer_survey_data, field_data, by="checklist_id")
# ebird_and_field_data$observation_count <- as.numeric(ebird_and_field_data$observation_count)
# 
# 
# 
# addFieldData <- function(bird_record) 
# {
#   print(bird_record)
#   print(as.data.frame(bird_record))
#   return(merge(, field_data %>% filter(point_count_checklist_id == as.character(bird_record["checklist_id"]))))
# }
# 
# ebird_and_field_data <- as.data.frame(apply(summer_survey_data, FUN=addFieldData, MARGIN=2))
# 
# 
# 
# 
# # NOTE need to add in code to load all_summaries and add wetness from field data etc from other script
# #  This presupposes those are already in memory, fix this for reproducibility
# audio_summary_data <- (arrange(all_summaries %>% group_by(wet, common_name) %>% summarize(best_confidence_median_audio = median(best_confidence), total_count_audio = sum(count), median_count_audio = median(count)), common_name, wet))
# point_count_summary_data <- (ebird_and_field_data %>% group_by(common_name, wet) %>% summarize(count = n(), total = sum(observation_count), median = median(observation_count)))
# 
# audio_and_point_counts <- merge(point_count_summary_data, audio_summary_data, by=c("wet","common_name"))
# 
# # species in both wet and dry environments
# wet_and_dry_species <- (audio_and_point_counts %>% group_by(common_name) %>% tally() %>% filter(n > 1))$common_name
# 
# wet_counts <- audio_and_point_counts %>% filter(wet==TRUE, common_name %in% wet_and_dry_species)
# dry_counts <- audio_and_point_counts %>% filter(wet==FALSE, common_name %in% wet_and_dry_species)
# wet_dry_difference <- wet_counts
# wet_dry_difference$best_confidence_median_audio <- wet_dry_difference$best_confidence_median_audio / dry_counts$best_confidence_median_audio
# wet_dry_difference$count <- wet_dry_difference$count / dry_counts$count
# wet_dry_difference$total <- wet_dry_difference$total / dry_counts$total
# wet_dry_difference$median <- wet_dry_difference$median / dry_counts$median
# wet_dry_difference$total_count_audio <- wet_dry_difference$total_count_audio / dry_counts$total_count_audio
# wet_dry_difference$best_confidence_median_audio <- wet_dry_difference$best_confidence_median_audio / dry_counts$best_confidence_median_audio
# wet_dry_difference$median_count_audio <- wet_dry_difference$median_count_audio / dry_counts$median_count_audio
# 
# 
# all_summaries_metadata <- merge(all_summaries, relevant_metadata[,c("checklist_string","point_count_checklist_id")] %>% mutate(checklist = checklist_string), by="checklist") %>%
#   mutate(checklist_id = point_count_checklist_id)
# 
# 
# 
# 
# # Get list of bird species on eBird from Mission / Rattlesnake Area
# ebird_mc_rc <- f_in %>% 
#   # 1. reference file
#   auk_ebd() %>% 
#   # 2. define filters
#   #  auk_species(species = "Oak Titmouse") %>% 
#   # 3. run filtering
#   auk_filter(file = f_out, overwrite=TRUE) %>% 
#   # 4. read text file into r data frame
#   read_ebd() %>%
#   filter(locality %in% c("Tunnel Trail","Rattlesnake Canyon Trail","Skofield Park","Santa Barbara Botanic Garden","Rocky Nook Park"))



