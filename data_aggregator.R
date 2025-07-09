
library(tidyverse)
library(janitor)
library(sf)

# Superdirectory containing directories for each checklist
#data_directory <- "D:/birdrec/data/birdrec_2022/madera_canyon/April_2022/"
#data_directory <- "D:/birdrec/data/birdrec_2022/sedgwick_ranch/"
#data_directory <- "D:/birdrec/data/goleta_summer_2022/san_antonio_creek/"
#data_directory <- "E:/santa_clara/bird_detections/test_data_02_16_2023/"
#data_directory <- "D:/birdrec/bird_detections/all_combined/"
#data_directory <- "G:/Bioacoustics/SCR_2023/detections/"
#compiled_directory <- "G:/Bioacoustics/SCR_2023/detections_compiled/"
data_directory <- "G:/Bioacoustics/Goleta_Various/detections/"
compiled_directory <- "G:/Bioacoustics/Goleta_Various/detections_compiled/"

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
  checklist_directory <- paste(data_directory,checklist,sep="")
  files <- list.files(path=checklist_directory)
  csv_files <- files[grepl(".*\\.csv", files)] 
  print(paste("Working on data from ", checklist, ", which has ", length(files), " detection files."), sep="")

  dataframe_list <- lapply(csv_files, ingestRecording, directory=checklist_directory, checklist_name=checklist)
  # Drop any empty dataframes (no birds detected)
  dataframe_list <- dataframe_list[lapply(dataframe_list, nrow) > 0]
  # Convert list of dataframes to one dataframe
  output_df <- bind_rows(dataframe_list)
  write_csv(output_df, paste(compiled_directory,"/",checklist,".csv",sep=""))
  return(output_df)
}
# Get all checklist IDs within data_directory
checklists <- list.files(path=data_directory)
# Drop any .csv readme or config files
checklists <- checklists[(1:length(checklists))*(1-grepl(".*\\.csv", checklists))]
# Collect all data from files
all_detections <- lapply(checklists, getBirdData)
# Remove any empty checklists
all_detections <- bind_rows(all_detections)
write_csv(all_detections, paste(compiled_directory, "all_detections.csv", sep=""))

# Summarize data - best detection confidence for each site
all_summaries <- all_detections %>% 
  group_by(site, common_name) %>% 
  summarize(best_confidence=max(confidence), 
            count=n()) %>%
  arrange(-best_confidence) %>%
  filter(best_confidence > 0.5)

view(arrange(all_summaries, common_name, site))

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

# Print out some data
write_csv(all_summaries, "C:/Users/grad/Downloads/all_summaries.csv")
write_csv(species_presence, "C:/Users/grad/Downloads/species_presence.csv")

# Generate diurnal signal for a bird
#   This automatically strips timestamps from filenames for SongMeter Micro and Audiomoth filenames
#   Assumes input is a list of dataframes, with separate dataframes for each site
getTimestampData <- function(detection_data)
{
  filenames <- detection_data$filename
  # IF the file is a Song Meter file
  if(TRUE)#length(grep("SMM", detection_data[1,]$filename)) > 0)
  {
    timestamps <- substr(filenames, nchar(filenames)-25, nchar(filenames)-20)
    dates <- substr(filenames, nchar(filenames)-34, nchar(filenames)-27)
  }
  # IF the file is NOT a SongMeter file
  else
  {
    timestamps <- substr(filenames, nchar(filenames)-25, nchar(filenames)-20)
    dates <- substr(filenames, nchar(filenames)-34, nchar(filenames)-27)
  }
  detection_data_mod <- detection_data
  # Set up time
  detection_data_mod$timestamp <- timestamps
  detection_data_mod$hour <- as.numeric(substr(detection_data_mod$timestamp, 1,2))
  detection_data_mod$minute <- as.numeric(substr(detection_data_mod$timestamp, 3,4))
  detection_data_mod$second <- as.numeric(substr(detection_data_mod$timestamp, 5,6))
  detection_data_mod$hour_fractional <- detection_data_mod$hour + detection_data_mod$minute/60 + detection_data_mod$second / 3600
  # Set up date
  detection_data_mod$date <- dates
  detection_data_mod$year <- as.numeric(substr(detection_data_mod$date, 1,4))
  detection_data_mod$month <- as.numeric(substr(detection_data_mod$date, 5,6))
  detection_data_mod$day <- as.numeric(substr(detection_data_mod$date, 7,8))
  detection_data_mod$yday <- lubridate::yday(paste(detection_data_mod$year,
                                                   detection_data_mod$month,
                                                   detection_data_mod$day,
                                                   sep="-"))
  
  return(detection_data_mod)
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
getSiteName <- function(datarow){return(relevant_metadata[relevant_metadata$checklist_string == datarow$checklist,]$site_name)}
all_summaries$wet <- rep(0,nrow(all_summaries))
all_summaries$site_name <- rep("",nrow(all_summaries))
for(ind in 1:nrow(all_summaries))
{
  all_summaries[ind,]$wet <- getWetness(all_summaries[ind,])
  all_summaries[ind,]$site_name <- as.character(getSiteName(all_summaries[ind,]))
}
all_summaries$creek <- substr(all_summaries$site_name, 1, nchar(all_summaries$site_name)-3)
all_summaries_mc_rc <- all_summaries %>% filter(creek %in% c("MC","RC"))

all_detections_timestamps$wet <- rep(NA,nrow(all_detections_timestamps))
all_detections_timestamps$checklist <- all_detections_timestamps$site
for(ind in 1:nrow(all_detections_timestamps))
{
  all_detections_timestamps[ind,]$wet <- getWetness(all_detections_timestamps[ind,])
}


# Get total volume of data collected at each point (which is proportional to total recording time):
sum(file.info(list.files(data_directory, all.files = TRUE, recursive = TRUE))$size)



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



