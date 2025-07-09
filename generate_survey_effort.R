
# Generate standardized .csv files which contain amount of survey effort at each site and time

# Load libraries
library(tidyverse)
library(here)
library(tuneR)

# Target directory which contains all the audio files
target_directory <- "G:/Bioacoustics/Goleta_2024/source/"
# Output directory to save effort files
effort_output_directory <- "G:/Bioacoustics/Goleta_2024/effort/"

generateSurveyEffort <- function(audio_filename)
{
  # Get date and time of audio file
  filename_chunks <- str_split(audio_filename, "/")[[1]]
  filename_without_path <- filename_chunks[[length(filename_chunks)]]
  datetime_str <- substring(filename_without_path, 
                            nchar(filename_without_path)-18, 
                            nchar(filename_without_path)-4)
  
  # Get Sitename
  sitename <- filename_chunks[[length(filename_chunks)-1]]
  
  # Set up metadata 
  metadata <- data.frame(site = sitename,
                         year = as.numeric(substr(datetime_str, 1,4)),
                         month = as.numeric(substr(datetime_str, 5,6)),
                         day = as.numeric(substr(datetime_str, 7,8)),
                         hour = as.numeric(substr(datetime_str, 10,11)),
                         minute = as.numeric(substr(datetime_str, 12,13)),
                         second = as.numeric(substr(datetime_str, 14,15))) %>% 
    mutate(doy = lubridate::yday(paste(year,month,day,sep="-")),
           time_frac = hour/24 + minute/24/60 + second/24/60/60)
  
  # Get the length of the audio file
  audioHeader <- readWave(audio_filename, header=TRUE)
  metadata$length <- round(audioHeader$samples / audioHeader$sample.rate, 2)
  
  # Write metadata to a file to record which minutes were recorded at each site
  survey_time_file <- paste(effort_output_directory, "survey_time.csv", sep="")
  write_csv(metadata, survey_time_file,
            append = file.exists(survey_time_file))
  
  return(metadata)
}



# Iterate across sites
site_dirs <- list.files(target_directory)
all_effort <- lapply(paste(target_directory, site_dirs, sep="/"), 
                         function(site_dir){
                           print(site_dir)
                           if(!dir.exists(site_dir))
                             return(-1)
                           filenames <- list.files(site_dir, pattern=".wav")
                           # Iterate across files within a site 
                           temp <- lapply(paste(site_dir, filenames, sep="/"),
                                          function(filename){
                                            return(generateSurveyEffort(filename))
                                          })
                           return(1)
                         })
