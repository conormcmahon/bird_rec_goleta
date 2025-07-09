
# For each year
#  For each site
#   For each 1 min audio file
#    For each detection
#     Get datetime
#     Get species numeric index (from allowed species list)
#     Get embedding 
#     Get site 
#     Append embedding to big fucking text file for that species, combined with an index
#     Append datetime, species index, site, and embedding ID into dataframe
# For each species
#  Open entire embeddings file in memory (prohibitively large?)
#  k-means clustering --> identify ~10 important classes (tune this number empirically per species)
#  human classifies clusters to important call types (e.g. song, chip, tsip, begging...)
#  For each 1 min audio file
#   generate a histogram of calls by type and confidence (include noise level and mockingbird, similar spp, human sound?)
#  Pick a subset of ~50 audio files for each combination of values 

library(tidyverse)
library(lubridate)

embeddings_dir <- "G:/Bioacoustics/Goleta_2023/embeddings"
detections_dir <- "G:/Bioacoustics/Goleta_2023/detections"
audio_chunk_spacing <- 2.5

all_detections <- rbind(read_csv("G:/Bioacoustics/Goleta_2023/detections/all_detections_compiled_2023.csv"),
                        read_csv("G:/Bioacoustics/Goleta_2022/detections/all_detections_compiled.csv"))

postprocessDetections  <- function(detections_file)
{
  # Get date and time of audio file
  filename_chunks <- str_split(detections_file, "/")[[1]]
  filename_without_path <- filename_chunks[[length(filename_chunks)]]
  datetime_str <- substring(filename_without_path, 1, 15)
  
  # Get Sitename
  sitename <- filename_chunks[[length(filename_chunks)-1]]
  
  # Generate site and datetime metadata
  metadata <- data.frame(site = sitename,
                         year = as.numeric(substr(datetime_str, 1,4)),
                         month = as.numeric(substr(datetime_str, 5,6)),
                         day = as.numeric(substr(datetime_str, 7,8)),
                         hour = as.numeric(substr(datetime_str, 10,11)),
                         minute = as.numeric(substr(datetime_str, 12,13)),
                         second = as.numeric(substr(datetime_str, 14,15))) %>% 
    mutate(doy = lubridate::yday(paste(year,month,day,sep="-")),
           time_frac = hour/24 + minute/24/60 + second/24/60/60)
  
  # Load detections table 
  detections <- read_csv(detections_file, show_col_types=FALSE) %>%
    janitor::clean_names()
  if(nrow(detections) == 0)
     return()
  
  # Add metadata to detections table
  detections_with_metadata <- merge(metadata, detections) %>% 
    mutate(time_frac = time_frac + start_s/24/60/60,
           audio_chunk = start_s/audio_chunk_spacing+1) %>% 
    select(-doy, -scientific_name, -end_s)
  
  return(detections_with_metadata) 
}

# Iterate across sites
site_dirs <- list.files(detections_dir)
all_detections <- lapply(paste(detections_dir, site_dirs, sep="/"), 
                         function(site_dir){
                           print(site_dir)
                           filenames <- list.files(site_dir)
                           # Iterate across files within a site 
                           newdata <- lapply(paste(site_dir, filenames, sep="/"),
                                             function(filename){
                                               return(postprocessDetections(filename))
                                             })
                           sitename <- str_split(site_dir, "/")[[1]]
                           sitename <- sitename[[length(sitename)]]
                           newdata_empty_df <- which(unlist(lapply(newdata, nrow)) == 0)
                           print(paste("  This folder contained ", length(newdata), " files.", sep=""))
                           if(length(newdata_empty_df > 0))
                             newdata <- bind_rows(newdata[-newdata_empty_df])
                           else
                             newdata <- bind_rows(newdata)
                           print(paste("  The total number of detections was ", nrow(newdata), sep=""))
                           data_summary <- (newdata %>% group_by(common_name) %>% summarize(count=n()) %>% arrange(-count))
                           print(paste("  The most common species was ",
                                       as.character(data_summary$common_name)[[1]],
                                       " with ", as.numeric(data_summary$count)[[1]], " detections.", sep=""))
                           write_csv(newdata, paste(detections_dir, "/", sitename, ".csv", sep=""))
                           return(newdata)
                         })
all_detections <- bind_rows(bind_rows(all_detections))
write_csv(all_detections, 
          paste(detections_dir, "all_detections_compiled.csv", sep="/"))

# Aggregate embeddings by species 


# Function which, for a given detection, gets the appropriate embedding, adds metadata, and 
#   prints to a csv which aggregates all embeddings for each species
associateEmbeddingAndMetadata <- function(ind, current_detections, current_embeddings)
{
  # Get time/location metadata for this detection
  target_detection <- current_detections[ind,]
  # Get embedding for this detection
  target_data <- current_embeddings[[target_detection$audio_chunk]]
  # Set up output dataframe combining metadata and embedding
  output_data <- data.frame(site = as.character(target_detection$site),
                            year = as.numeric(target_detection$year),
                            month = as.numeric(target_detection$month),
                            day = as.numeric(target_detection$day),
                            hour = as.numeric(target_detection$hour),
                            minute = as.numeric(target_detection$minute),
                            second = as.numeric(target_detection$second),
                            start_s = as.numeric(target_detection$start_s),
                            time_frac = as.numeric(target_detection$time_frac),
                            confidence = as.numeric(target_detection$confidence))
  num_metadata_columns <- ncol(output_data)
  #  Add embedding data and correct the names
  output_data[1,(num_metadata_columns+1):(num_metadata_columns+1024)] <- target_data
  names(output_data) <- c(names(output_data)[1:num_metadata_columns], sprintf("emb_%04d", 1:1024))
  # Check whether output directory exists, and create it if not
  output_dir <- paste(embeddings_dir, 
                      "..",
                      "species_embeddings", 
                      sep="/")
  dir.create(file.path(output_dir), showWarnings = FALSE)
  # Set up filepath for output, based on species name
  output_file <- paste(output_dir, "/",
                       target_detection$common_name %>% janitor::make_clean_names(),
                       ".csv", sep="")
  # Write file 
  write_csv(output_data, 
            output_file,
            append = file.exists(output_file))
  
  return(output_data)
}

postprocessEmbeddings <- function(embeddings_filename)
{
  # Get date and time of audio file
  filename_chunks <- str_split(embeddings_filename, "/")[[1]]
  filename_without_path <- filename_chunks[[length(filename_chunks)]]
  datetime_str <- substring(filename_without_path, 1, 15)
  
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
  
  # Write metadata to a file to record which minutes were recorded at each site
  survey_time_file <- paste(detections_dir, "/../survey_time.csv", sep="")
  write_csv(metadata, survey_time_file,
            append = file.exists(survey_time_file))
  
  # Check whether there were any detections in this segment
  detections <- all_detections %>% 
    filter(site == sitename,
           year == metadata$year, 
           month == metadata$month, 
           day == metadata$day, 
           hour == metadata$hour, 
           minute == metadata$minute, 
           second == metadata$second)
  # If there were no detections, skip this embeddings file
  if(nrow(detections) == 0)
    return()
  
  # Otherwise, load embedding data for this audio file 
  # Load text file and get vector in string CSV format
  embeddings_text <- read.delim(embeddings_filename, header=FALSE)[,3]
  # Embeddings are split into 3-second chunks; gather them all into a list of embedding vectors 
  embeddings <- lapply(embeddings_text, 
                       function(data_line) {
                         as.numeric(strsplit(data_line, ",")[[1]])
                       })
  rm(embeddings_text)
  
  return(lapply(1:nrow(detections),
                associateEmbeddingAndMetadata,
                current_detections=detections,
                current_embeddings=embeddings))
}

# Iterate across sites
site_dirs <- list.files(embeddings_dir)
all_embeddings <- lapply(paste(embeddings_dir, site_dirs, sep="/"), 
                         function(site_dir){
                           print(site_dir)
                           if(!dir.exists(site_dir))
                             return(-1)
                           filenames <- list.files(site_dir)
                           # Iterate across files within a site 
                           temp <- lapply(paste(site_dir, filenames, sep="/"),
                                                    function(filename){
                                                      return(postprocessEmbeddings(filename))
                                                    })
                           # newdata <- lapply(paste(site_dir, filenames, sep="/"),
                           #                   function(filename){
                           #                     return(postprocessEmbeddings(filename))
                           #                   })
                           # newdata_empty_df <- which(unlist(lapply(newdata, nrow)) == 0)
                           # print(paste("  This folder contained ", length(newdata), " embeddings files.", sep=""))
                           # if(length(newdata_empty_df > 0))
                           #   newdata <- bind_rows(newdata[-newdata_empty_df])
                           # else
                           #   newdata <- bind_rows(newdata)
                           return(1)
                         })

