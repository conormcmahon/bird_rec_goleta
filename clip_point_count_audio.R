
# Load required packages
library(readr)
library(tuneR)
library(sf)
library(stringr)

# Pseudocode
# 
# Iterate over all points in the field data
# Skip points with no eBird checklist for the point count
# For points which did have a point count, get the timestamp (day hour minute second)
# Get the SD card name and site name 
# Get the filepath for the audio file
# Check that the audio file exists... if not, flag and skip this point
# Determine the second of the audio file to begin the sample based on filename and field timestamp
# Bump back 1 minute, and clip 7 minutes (1 minute buffer on each side to account for drift) 
# Load the audio file
# Write the output somewhere as .wav, with the site name and date in the filename

# Read field data file
field_data <- st_read(here::here("spatial_data","birdrec_site_goleta_2022.gpkg"))

# For some reason, some of the timestamps are being read incorrectly out of the gpkg
#   It seems like this is a time zone issue. All the incorrect stamps appear to be 7 hours early
#   Our timezone is UTC +7, so it must be related to this
#   Only some timestamps are affected and I'm not sure why
#   For now, just manually correcting this in a hacky way:
point_count_timestamp_corrected <- field_data$point_count_timestamp
broken_timestamps <- which((as.numeric(format(point_count_timestamp_corrected, "%H")) < 6))
point_count_timestamp_corrected[broken_timestamps] <- point_count_timestamp_corrected[broken_timestamps] + 3600*7
field_data$point_count_timestamp <- point_count_timestamp_corrected


# Set up mapping between creek abbreviation code and longform creek name
#    long form names
site_names <- c("arroyo_quemado",
                "atascadero_creek",
                "ellwood_grove",
                "kinevan_road",
                "lake_los_carneros",
                "mission_creek",
                "atascadero_creek",
                "ncos",
                "rattlesnake_creek",
                "san_antonio_creek",
                "san_jose_creek",
                "sedgwick")
#    site name codes
names(site_names) <- c("BR",
                       "AC",
                       "EW",
                       "KR",
                       "LLC",
                       "MC",
                       "MYC",
                       "NCOS",
                       "RC",
                       "SAC",
                       "SJC",
                       "SR")


# Generate diurnal signal for a bird
#   This automatically strips timestamps from filenames for SongMeter Micro and Audiomoth filenames
#   Assumes input is a list of dataframes, with separate dataframes for each site
getFileTimestamp <- function(filepath)
{
  filename <- str_split(filepath, "/")
  filename <- filename[[1]][length(filename[[1]])]
  # IF the file is a Song Meter file
  if(length(grep("SMM", filename)) > 0)
  {
    timestamp_string <- substr(filename, nchar(filename)-18, nchar(filename)-4)
  }
  # IF the file is NOT a SongMeter file
  else
  {
    timestamp_string <- substr(filename, nchar(filename)-18, nchar(filename)-4)
  }
  
  # Build a POSIXct version of timestamp
  # format:  "YYYY-MM-DD HH:MM:SS PDT"
  file_timestamp <- paste(substr(timestamp_string, 1,4), "-",
                          substr(timestamp_string, 5,6), "-",
                          substr(timestamp_string, 7,8), " ",
                          substr(timestamp_string, 10,11), ":",
                          substr(timestamp_string, 12,13), ":",
                          substr(timestamp_string, 14,15), " PDT",
                          sep="")
  file_timestamp <- as.POSIXct(file_timestamp)
  
  return(file_timestamp)
}


# For a given row from the field dataframe, find the corresponding audio file and subset the portion from the checklist
subsetPointCountAudio <- function(data_row, checklist_length, time_buffer, output_directory)
{
  # SD card numeric ID
  sd_card <- as.character(data_row$sd_card_id)
  # creek abbreviation from site code (e.g. LLC for Lake Los Carneros)
  site_abbreviation <- as.character(unlist(lapply(str_split(unique(data_row$site_name), '-'), function(x){x[[1]]})))
  # lookup longform site name in map
  site_name <- as.character(site_names[site_abbreviation])
  # skip Kinevan Road, because we didn't record there
  if(site_name == "kinevan_road")
    return(-1)
  if(is.na(data_row$point_count_timestamp))
    return(-1)
  # Date of unit collection (used in sorting audio data)
  collection_year <- format(data_row$timestamp_collected, '%Y')
  collection_month <- format(data_row$timestamp_collected, '%m')
  collection_day <- format(data_row$timestamp_collected, '%d')
  
  # Path to folder containing audio files for this site
  audio_folder_path <- here::here("..","..","data","goleta_summer_2022",site_name,paste(paste("sd",sd_card,collection_year,collection_month,collection_day,sep="_"),"/",sep=""))
  
  # List of all audio files in that folder  
  files <- list.files(audio_folder_path)
  
  # Get timestamps out of all filenames
  file_timestamps <- (lapply(files, getFileTimestamp))
  
  # Compare point count time to audio file timestamps
  time_difference <- unlist(as.numeric(lapply(file_timestamps, function(x){difftime(data_row$point_count_timestamp, x, units="secs")})))
  
  # If the point count was before all timestamps, then give up
  if(sum(time_difference < 0) == length(time_difference))
  {
    print(paste("For site ", data_row$site_name, " and SD card ", sd_card, ", the point count was before recording started."), sep="")
    print(time_difference)
    return(-2)
  }
  
  # Find the last audio recording which started before the point count
  post_count_file <- NA
  post_count_time_difference <- NA
  if(sum(time_difference > 0) == length(time_difference))
  {
    pre_count_file <- files[length(files)]
    pre_count_time_difference <- time_difference[length(time_difference)]
  }
  else
  {
    # Which audio file timestamps are before the point counts?
    all_pre_count_files <- which(time_difference > 0)
    # And the FIRST audio file which starts before the count
    pre_count_file <- files[all_pre_count_files[length(all_pre_count_files)]]
    pre_count_time_difference <- time_difference[all_pre_count_files[length(all_pre_count_files)]]
    # And the NEXT audio file - the first which started after the audio file
    post_count_file <- files[all_pre_count_files[length(all_pre_count_files)]+1]
    post_count_time_difference <- time_difference[all_pre_count_files[length(all_pre_count_files)]+1]
  }
  print(pre_count_file)
  print(data_row$point_count_timestamp)
  print(pre_count_time_difference)
  # Read in the audio file during which the checklist started
  pre_count_audio <- readWave(paste(audio_folder_path, pre_count_file, sep="/"))
  print(pre_count_audio)
  # Record the length (in seconds) of that file
  pre_count_audio_length <- round(length(pre_count_audio@left) / pre_count_audio@samp.rate, 2)
  print(pre_count_audio_length)
  # Check whether the file is long enough
  #   If it's not, that's a sign that the recorder stopped before the point count started (ran out of battery, storage, etc.)
  if(pre_count_audio_length < pre_count_time_difference)
  {
    print(paste("For site ", data_row$site_name, " and SD card ", sd_card, ", the point count started after the recorder stopped recording."), sep="")
    return(-3)
  }
  # Record the time stamp (in minutes) at which the 
  checklist_start_time <- pre_count_time_difference - time_buffer
  checklist_end_time <- checklist_start_time + checklist_length + time_buffer*2
  post_count_checklist_end_time <- NA
  if(checklist_end_time > pre_count_audio_length)
  {
    if(!is.na(post_count_file))
    {
      post_count_checklist_end_time <- checklist_end_time - pre_count_audio_length
    }
    checklist_end_time <- pre_count_audio_length
  }
  print(checklist_start_time)
  print(checklist_end_time)
  
  subset_audio <- readWave(paste(audio_folder_path, pre_count_file, sep="/"), from=checklist_start_time, to=checklist_end_time, units="seconds")
  if(!is.na(post_count_checklist_end_time))
  {
    subset_audio_backend <- readWave(paste(audio_folder_path, post_count_file, sep="/"), from=0, to=post_count_checklist_end_time, units="seconds")
    subset_audio <- bind(subset_audio, subset_audio_backend)
  }

  writeWave(subset_audio, paste(output_directory,"/",data_row$site_name,"_",collection_year,collection_month,collection_day,".wav",sep=""))
  
  return(1)
}


output_codes <- rep(0, nrow(field_data))
for(ind in 1:nrow(field_data))
{
  output_codes[ind] <- subsetPointCountAudio(field_data[ind,], checklist_length=300, time_buffer=60, output_directory=here::here("..","..","data","point_count_audio"))
}
field_data$output_code <- output_codes