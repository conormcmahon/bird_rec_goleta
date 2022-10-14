
library(tidyverse)
library(janitor)
library(sf)

# Superdirectory containing directories for each checklist
#data_directory <- "D:/birdrec/data/birdrec_2022/madera_canyon/April_2022/"
#data_directory <- "D:/birdrec/data/birdrec_2022/sedgwick_ranch/"
#data_directory <- "D:/birdrec/data/goleta_summer_2022/sedgwick/"
#data_directory <- "D:/birdrec/data/goleta_summer_2022/san_antonio_creek/"
data_directory <- "D:/birdrec/bird_detections/all_combined/"

# Ingest a single recording from BirdNet .csv tab-delimited report
ingestRecording <- function(filename, directory, checklist_name)
{
  data <- read_csv(paste(directory,filename,sep="/"),
                   col_types = cols(start = col_character(),
                                    end = col_character(),
                                    scientific_name = col_character(),
                                    common_name = col_character(),
                                    confidence = col_character())) %>% 
    janitor::clean_names()
  data$filename <- rep(paste(directory,filename,sep="/"), nrow(data))
  data$site <- rep(checklist_name, nrow(data))
  return(data)
}
# Process all recording files associated with a single checklist directory
getBirdData <- function(checklist)
{
  checklist_directory <- paste(data_directory,checklist,sep="/")
  files <- list.files(path=checklist_directory)
  csv_files <- files[grepl(".*\\.csv", files)] 

  dataframe_list <- lapply(csv_files, ingestRecording, directory=checklist_directory, checklist_name=checklist)
  # Drop any empty dataframes (no birds detected)
  dataframe_list <- dataframe_list[lapply(dataframe_list, nrow) > 0]
  # Convert list of dataframes to one dataframe
  bind_rows(dataframe_list)
}
# Get all checklist IDs within data_directory
checklists <- list.files(path=data_directory)
# Drop any .txt readme or config files
checklists <- checklists[(1:length(checklists))*(1-grepl(".*\\.txt", checklists))]
# Collect all data from files
all_detections <- lapply(checklists, getBirdData)
# Remove any empty checklists
all_detections <- all_detections[(unlist(lapply(all_detections,nrow))>0)] 


# Summarize data - best detection confidence for each checklist
#   Function to apply for each checklist
bestDetections <- function(detections)
{
  best_detection <- detections %>% 
    group_by(common_name) %>% 
    summarize(best_confidence=max(confidence), 
              count=n()) %>%
    arrange(-best_confidence) %>%
    filter(best_confidence > 0.5)
} 
#   Apply the above to each checklist
checklist_summaries <- lapply(all_detections, bestDetections)
names(checklist_summaries) <- checklists[1:length(checklist_summaries)]

all_summaries <- bind_rows(checklist_summaries, .id="checklist")
view(all_summaries)

# Get a list of all species detected with high confidence, and also
#   total - total number of detections across all checklists
#   presence_count - number of checklists on which the species was detected
#   presence_rate - 
species_presence <- all_summaries %>% 
  group_by(common_name) %>% 
  summarize(total=sum(count),
            presence_absolute=n(),
            presence_rate=n()/length(checklists),
            best_confidence=max(best_confidence)) %>% 
  arrange(-total)
view(species_presence)

# Generate diurnal signal for a bird
#   This automatically strips timestamps from filenames for SongMeter Micro and Audiomoth filenames
#   Assumes input is a list of dataframes, with separate dataframes for each site
getTimestampData <- function(detection_data_list)
{
  getTimeStamps <- function(detection_data)
  {
    filenames <- detection_data$filename
    # IF the file is a Song Meter file
    if(length(grep("SMM", detection_data[1,]$filename)) > 0)
    {
      timestamps <- substr(filenames, nchar(filenames)-25, nchar(filenames)-20)
    }
    # IF the file is NOT a SongMeter file
    else
    {
      timestamps <- substr(filenames, nchar(filenames)-25, nchar(filenames)-20)
    }
    detection_data_mod <- detection_data
    detection_data_mod$timestamp <- timestamps
    detection_data_mod$hour <- as.numeric(substr(detection_data_mod$timestamp, 1,2))
    detection_data_mod$minute <- as.numeric(substr(detection_data_mod$timestamp, 3,4))
    detection_data_mod$second <- as.numeric(substr(detection_data_mod$timestamp, 5,6))
    detection_data_mod$hour_fractional <- detection_data_mod$hour + detection_data_mod$minute/60 + detection_data_mod$second / 3600
    
    return(detection_data_mod)
  }
  
  data_list_with_stamps <- lapply(detection_data_list, getTimeStamps)
  
  data_out <- bind_rows(data_list_with_stamps)
}

all_detections_timestamps <- getTimestampData(all_detections)

all_detections_timestamps

species = "California Scrub-Jay"
hist((all_detections_timestamps %>% filter(common_name == species, confidence>0.1))$hour_fractional, xlab="Hour of Day", ylab="Number of Calls", main=paste("Daily Vocal Activity - ",species,sep=""), breaks=0:24)
print(all_summaries %>% filter(common_name == species))
print(species_presence %>% filter(common_name == species))

species = "House Finch"
hist((all_detections_timestamps %>% filter(common_name == species, confidence>0.1))$hour_fractional, xlab="Hour of Day", ylab="Number of Calls", main=paste("Daily Vocal Activity - ",species,sep=""), breaks=0:24)
print(all_summaries %>% filter(common_name == species))
print(species_presence %>% filter(common_name == species))

species = "California Quail"
hist((all_detections_timestamps %>% filter(common_name == species, confidence>0.1))$hour_fractional, xlab="Hour of Day", ylab="Number of Calls", main=paste("Daily Vocal Activity - ",species,sep=""), breaks=0:24)
print(all_summaries %>% filter(common_name == species))
print(species_presence %>% filter(common_name == species))

# Filter to species with relatively high confidence
high_confidence_species_presence <- species_presence %>% 
  filter(total > 150,
         best_confidence > 0.8)

# Load spatial data
spatial_data <- st_read(here::here("spatial_data","birdrec_site_goleta_2022.gpkg"))
# Clean up SD card names (force to be two digits, make character, remove SD and other noise)
# sd_card_ids <- as.character(
#                               unlist(lapply(str_extract_all(as.character(spatial_data$sd_card_id), "\\d*"), 
#                                             max)))
# forceTwoDigits <- function(sd_name)
# {
#   if(nchar(sd_name) == 1)
#   {
#     return(paste("0",sd_name,sep=""))
#   }
#   else if(nchar(sd_name == 2))
#   {
#     return(sd_name)        
#   }
#   else
#   {
#     print(paste("ERROR - this SD card name doesn't make sense: ",sd_name, sep=""))
#     return(sd_name)
#   }
# }
# sd_card_ids <- unlist(lapply(sd_card_ids, forceTwoDigits))
# spatial_data$sd_cleaned <- sd_card_ids

# Create Lookup from checklists to sites 
checklist_dates <- as.Date(paste(substr(checklists,7,10),substr(checklists,12,13),substr(checklists,15,16),sep="-"))
sd_card_ids <- substr(checklists,4,5)
getSpatialMatch <- function(ind)
{
  site_matches <- (as.character(spatial_data$sd_card_id) == sd_card_ids[ind]) * (as.Date(spatial_data$timestamp_collected) == checklist_dates[ind])
  site_matches[is.na(site_matches)] <- 0
  return (which(as.logical(site_matches)))
}

# replace NA in surface water values for spatial_data with 0
spatial_data[which(is.na(spatial_data$surface_flow)),]$surface_flow <- FALSE
spatial_data[which(is.na(spatial_data$pooling_water)),]$pooling_water <- FALSE
spatial_data[which(is.na(spatial_data$seeps_or_springs)),]$seeps_or_springs <- FALSE
# Add "wet" variable for presence of any surface water
spatial_data$wet <- (spatial_data$surface_flow + spatial_data$pooling_water + spatial_data$seeps_or_springs) > 0 

checklist_matches <- unlist(lapply(1:length(checklist_dates), getSpatialMatch))

# function to return a spatial_data row given a checklist string
getChecklistMetadata <- function(checklist_name)
{
  checklist_date <- as.Date(paste(substr(checklist_name,7,10),substr(checklist_name,12,13),substr(checklist_name,15,16),sep="-"))
  sd_card_id <- substr(checklist_name,4,5)
  
  site_matches <- (as.character(spatial_data$sd_card_id) == sd_card_id) * (as.Date(spatial_data$timestamp_collected) == checklist_date)
  site_matches[is.na(site_matches)] <- 0
  return (spatial_data[which(as.logical(site_matches)),])
}

relevant_metadata <- (bind_rows(lapply(checklists, getChecklistMetadata)))
relevant_metadata$checklist_string <- checklists

getWetness <- function(datarow){return(relevant_metadata[relevant_metadata$checklist_string == datarow$checklist,]$wet)}
all_summaries$wet <- rep(0,nrow(all_summaries))
for(ind in 1:nrow(all_summaries))
{
  all_summaries[ind,]$wet <- getWetness(all_summaries[ind,])
}

all_detections_timestamps$wet <- rep(NA,nrow(all_detections_timestamps))
all_detections_timestamps$checklist <- all_detections_timestamps$site
for(ind in 1:nrow(all_detections_timestamps))
{
  all_detections_timestamps[ind,]$wet <- getWetness(all_detections_timestamps[ind,])
}



yewa_sd_ids <- substr((all_summaries %>% filter(common_name == "Yellow Warbler"))$checklist, 4,5)
yewa_data <- (all_summaries %>% filter(common_name == "Yellow Warbler"))
yewa_data$geometry <- baron_sites[unlist(lapply(yewa_sd_ids, function(x) which(str_detect(x, baron_sites$sd_cleaned)))),]$geometry
st_write(yewa_data, here::here("spatial_data","yewa.gpkg"))

psfl_sd_ids <- substr((all_summaries %>% filter(common_name == "Pacific-slope Flycatcher"))$checklist, 4,5)
psfl_data <- (all_summaries %>% filter(common_name == "Pacific-slope Flycatcher"))
psfl_data$geometry <- baron_sites[unlist(lapply(psfl_sd_ids, function(x) which(str_detect(x, baron_sites$sd_cleaned)))),]$geometry
st_write(psfl_data, here::here("spatial_data","psfl.gpkg"))

wiwa_sd_ids <- substr((all_summaries %>% filter(common_name == "Wilson's Warbler"))$checklist, 4,5)
wiwa_data <- (all_summaries %>% filter(common_name == "Wilson's Warbler"))
wiwa_data$geometry <- baron_sites[unlist(lapply(wiwa_sd_ids, function(x) which(str_detect(x, baron_sites$sd_cleaned)))),]$geometry
st_write(wiwa_data, here::here("spatial_data","wiwa.gpkg"))

