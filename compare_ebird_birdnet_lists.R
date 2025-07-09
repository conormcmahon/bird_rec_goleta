

# This script compares the BirdNET Species list to an eBird EBD database file 
#   Birds which are listed in eBird but NOT the species list are noted
#   These can then be considered for inclusion to the BirdNET species list (optionally)

library(auk)
library(tidyverse)
library(here)
library(sf)

# path to the ebird data file, here all Santa Barbara County data up until July 2022
# get the path to the example data included in the package
f_in <- here::here("..","..","eBird","ebird_ca_dec_2023","ebd_US-CA_smp_relOct-2023.txt")
# output text file
f_out <- "all_sb_co_data_full_timeseries.txt"
ebird_data <- f_in %>% 
  # 1. reference file
  auk_ebd() %>% 
  # 2. define filters
  auk_county("US-CA-083") %>%
  # 3. run filtering
  auk_filter(file = f_out, overwrite=TRUE) %>% 
  # 4. read text file into r data frame
  read_ebd()
# Get summary of eBird species by number of total detections
ebird_species <- ebird_data %>% 
  group_by(common_name) %>%
  summarize(num_reports = n()) %>%
  arrange(-num_reports)


# Load BirdNET data
#   this comes in format    Genus species_Common Name
#   reformat to extract these 3 pieces of information separately
birdnet_df <- read_csv("D:/birdrec/protocol/riparian_broad_species_list/species_list.txt", col_names = "full_string") %>%
  mutate(binomial = unlist(lapply(str_split(full_string, "_"), function(newdata){return(newdata[[1]])})),
         common_name = unlist(lapply(str_split(full_string, "_"), function(newdata){return(newdata[[2]])}))) %>%
  mutate(genus = unlist(lapply(str_split(binomial, "_"), function(newdata){return(newdata[[1]])})),
         species = unlist(lapply(str_split(binomial, "_"), function(newdata){return(newdata[[1]])})))

write_csv(birdnet_df, "D:/birdrec/protocol/birdnet_riparian_species_list.xlsx")

# Look for species which are in eBird but not the BirdNET list
missing_from_birdnet <- ebird_species[!(ebird_species$common_name %in% birdnet_df$common_name),] %>%
  arrange(-num_reports)
View(missing_from_birdnet)

write_csv(missing_from_birdnet, "D:/birdrec/protocol/missing_from_birdnet_riparian_broad.csv")
