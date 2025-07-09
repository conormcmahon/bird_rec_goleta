
# Code to extract clusters of vocal types from feature embeddings of bird vocalizations identified to species 

library(tidyverse)  # handle dataframes
library(tuneR)      # load and subsample audio files to 3-sec snippets
library(sf)         # load field vector data (gpkg) to get unit type for each site

# Load field data (includes which unit type was used at each site)
field_data <- st_read("D:/birdrec/field_data_goleta_2022/goleta_qfield_final_2023/birdrec_site_goleta_2022.gpkg")
#  set up dictionary relating sitenames to ARU types 
aru_types <- substr(field_data$recording_unit_id, 1, 2)
names(aru_types) <- str_replace(field_data$site_name, "-", "_")

# Audio Directory
audio_directory <- "G:/Bioacoustics/Goleta_2023/source"

# Length in seconds of an audio snippet for one BirdNET detection
length_s <- 3

# Code to reconstruct filepaths to audiofiles based on metadata in dataframe
getAudioFilepath <- function(metadata, audio_dir, common_name)
{
  # Construct date and time for filename from metadata
  filename_bases <- paste(sprintf("%04d", metadata$year),
                         sprintf("%02d", metadata$month),
                         sprintf("%02d", metadata$day),
                         "_",
                         sprintf("%02d", metadata$hour),
                         sprintf("%02d", metadata$minute),
                         sprintf("%02d", metadata$second), 
                         sep="")
  # Double-check that year is correct in audio path
  audio_dir <- str_replace_all(audio_dir, "2022", as.character(metadata$year))
  audio_dir <- str_replace_all(audio_dir, "2023", as.character(metadata$year))
  # Append datetime to directory structure to reconstruct audio filepath
  original_paths <- paste(audio_dir, "/", metadata$site, "/", filename_bases, ".wav", sep="")
  # Create new filename to save test 3-sec audio snippet to 
  new_filenames <- paste(audio_dir, "/../call_segments/", common_name, "/", "cluster_", metadata$cluster, "/", metadata$site, "_", filename_bases, "_", metadata$start_s, sep="")
  
  return(data.frame(new_filename = new_filenames,
                    filename_base = filename_bases,
                    original_path = original_paths,
                    confidence = metadata$confidence,
                    cluster = metadata$cluster,
                    start_s = metadata$start_s,
                    unit_type = metadata$unit_type,
                    common_name = common_name))
}
# Copy 3 second audio snippet for a given bird detection to example folder
copyAudioSegment <- function(embeddings_subset_df)
{
  for(ind in 1:nrow(embeddings_subset_df))
  {
    # Read source file
    audio <- readWave(embeddings_subset_df[ind,"original_path"])
    
    # the frequency of your audio file
    freq <- audio@samp.rate
    # the length and duration of your audio file
    totlen <- length(audio)
    totsec <- totlen/freq 
    
    # Get target 3-second subsample
    subsamp <- matrix(audio@left, ncol=1, nrow=length(audio@left))
    start_i <- embeddings_subset_df[ind,"start_s"]*freq
    end_i <- min((embeddings_subset_df[ind,"start_s"]+length_s)*freq,
                 length(audio@left))
    subsamp <- subsamp[start_i:end_i]
    # Generate new Audio dataset
    audio_out <- Wave(left = subsamp, 
                      samp.rate = audio@samp.rate,
                      bit = audio@bit)
    # Check whether target directory exists, and create it if not
    slash_indices <- gregexpr("/", embeddings_subset_df[ind,"new_filename"])
    output_dir <- substr(embeddings_subset_df[ind,"new_filename"], 
                         1, slash_indices[[1]][length(slash_indices[[1]])])
    print(paste("  printing segment for ", output_dir, sep=""))
    if(!(dir.exists(output_dir)))
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    # Write to disk
    writeWave(audio_out, embeddings_subset_df[ind,"new_filename"])
  }
}

# Generate clusters and save example audio for a given species and recording unit type 
processSpecies <- function(common_name, num_clusters, unit, copy_audio_segment=TRUE)
{
  # Set up clean common name string for directory structure
  clean_name <- janitor::make_clean_names(common_name)
  
  # Load all embeddings for a target species 
  embeddings_all <- rbind(read_csv(paste("G:/Bioacoustics/Goleta_2023/species_embeddings/", clean_name, ".csv", sep="")),
                          read_csv(paste("G:/Bioacoustics/Goleta_2022/species_embeddings/", clean_name, ".csv", sep="")))
  #  assign ARU types to embeddings dataframe
  embeddings_all$unit_type <- aru_types[embeddings_all$site]
  
  # Pull out embeddings for target species and recording unit type
  target_embeddings <- embeddings_all %>% filter(unit_type == unit)
  # Check that there are enough embeddings in the data subset
  if(nrow(target_embeddings %>% filter(confidence > 0.5)) < 10*num_clusters)
  {
    print(paste("ERROR - not enough data in subset for ", common_name, " and ", unit, " recorders! Only ",
                nrow(target_embeddings %>% filter(confidence > 0.5)), " detections.", sep=""))
  }
  
  # Cluster on embeddings
  #  Set seed for RNG
  set.seed(1)
  #  Build X clusters, using ONLY the data with confidence > 0.5
  clusters_hc <- kmeans(target_embeddings %>% filter(confidence >= 0.5) %>% select(contains("emb")), centers=num_clusters)
  #  Assign clusters to all data, including even those with low confidence
  clusters_all <- kmeans(target_embeddings %>% select(contains("emb")), centers=clusters_hc$centers)
  target_embeddings$cluster <- clusters_all$cluster
  
  if(copy_audio_segment)
  {
    for(target_cluster in 1:num_clusters)
    {
      example_detections <- (getAudioFilepath(head(target_embeddings %>% 
                                                     select(-contains("emb")) %>% 
                                                     filter(cluster == target_cluster) %>% 
                                                     arrange(-confidence), 
                                                   100), 
                                              "G:/Bioacoustics/Goleta_2023/source",
                                              paste(clean_name, "/", unit, sep="")))
      copyAudioSegment(example_detections)
    }
  }
    
  return(target_embeddings)
}

# Generate clusters and save example snippets for AM and SM recordings of Yellow Warbler
yewa_am <- processSpecies("Yellow Warbler", 5, "AM", TRUE)
yewa_sm <- processSpecies("Yellow Warbler", 5, "SM", TRUE)

# Generate clusters and save example snippets for AM and SM recordings of Yellow Warbler
wiwa_am <- processSpecies("Wilson's Warbler", 5, "AM", TRUE)
wiwa_sm <- processSpecies("Wilson's Warbler", 5, "SM", TRUE)

# Generate clusters and save example snippets for AM and SM recordings of Yellow Warbler
wavi_am <- processSpecies("Warbling Vireo", 5, "AM", TRUE)
wavi_sm <- processSpecies("Warbling Vireo", 5, "SM", TRUE)

# Generate clusters and save example snippets for AM and SM recordings of Yellow Warbler
cbch_am <- processSpecies("Chestnut-backed Chickadee", 5, "AM", TRUE)
cbch_sm <- processSpecies("Chestnut-backed Chickadee", 5, "SM", TRUE)

# Generate clusters and save example snippets for AM and SM recordings of Yellow Warbler
pufi_am <- processSpecies("Purple Finch", 5, "AM", TRUE)
pufi_sm <- processSpecies("Purple Finch", 5, "SM", TRUE)

# Generate clusters and save example snippets for AM and SM recordings of Yellow Warbler
hawo_am <- processSpecies("Hairy Woodpecker", 5, "AM", TRUE)
hawo_sm <- processSpecies("Hairy Woodpecker", 5, "SM", TRUE)

# Generate clusters and save example snippets for AM and SM recordings of Yellow Warbler
coye_am <- processSpecies("Common Yellowthroat", 5, "AM", TRUE)
coye_sm <- processSpecies("Common Yellowthroat", 5, "SM", TRUE)

# Generate clusters and save example snippets for AM and SM recordings of Yellow Warbler
sosp_am <- processSpecies("Song Sparrow", 5, "AM", TRUE)
sosp_sm <- processSpecies("Song Sparrow", 5, "SM", TRUE)

# Generate clusters and save example snippets for AM and SM recordings of Acorn Woodpecker
acwo_am <- processSpecies("Acorn Woodpecker", 5, "AM", TRUE)
acwo_sm <- processSpecies("Acorn Woodpecker", 5, "SM", TRUE)

# Generate clusters and save example snippets for AM and SM recordings of Acorn Woodpecker
ocwa_am <- processSpecies("Orange-crowned Warbler", 5, "AM", TRUE)
ocwa_sm <- processSpecies("Orange-crowned Warbler", 5, "SM", TRUE)

# Generate clusters and save example snippets for AM and SM recordings of Acorn Woodpecker
spto_am <- processSpecies("Spotted Towhee", 5, "AM", TRUE)
spto_sm <- processSpecies("Spotted Towhee", 5, "SM", TRUE)



# Manually look at clusters --> determine which are calls vs. songs
# YEWA - AudioMoth
yewa_calltypes_am <- c("Song", "Song", "Call", "Song", "Song")
names(yewa_calltypes_am) <- 1:5
yewa_am$calltype <- as.character(yewa_calltypes_am[yewa_am$cluster])
# YEWA - Song Meter Micro
yewa_calltypes_sm <- c("Song", "Call", "Song", "Song", "Song")
names(yewa_calltypes_sm) <- 1:5
yewa_sm$calltype <- as.character(yewa_calltypes_sm[yewa_sm$cluster])

# WIWA - AudioMoth
wiwa_calltypes_am <- c("Song", "Song", "Call", "Song", "Song")
names(wiwa_calltypes_am) <- 1:5
wiwa_am$calltype <- as.character(wiwa_calltypes_am[wiwa_am$cluster])
# WIWA - Song Meter Micro
wiwa_calltypes_sm <- c("Song", "Song", "Call", "Song", "Song")
names(wiwa_calltypes_sm) <- 1:5
wiwa_sm$calltype <- as.character(wiwa_calltypes_sm[wiwa_sm$cluster])

# WAVI - AudioMoth
wavi_calltypes_am <- c("Song", "Song", "Call", "Song", "Song")
names(wavi_calltypes_am) <- 1:5
wavi_am$calltype <- as.character(wavi_calltypes_am[wavi_am$cluster])
# WAVI - Song Meter Micro
wavi_calltypes_sm <- c("Song", "Song", "Call", "Song", "Song")
names(wavi_calltypes_sm) <- 1:5
wavi_sm$calltype <- as.character(wavi_calltypes_sm[wavi_sm$cluster])


# WAVI - AudioMoth
wavi_calltypes_am <- c("Song", "Song", "Call", "Song", "Song")
names(wavi_calltypes_am) <- 1:5
wavi_am$calltype <- as.character(wavi_calltypes_am[wavi_am$cluster])
# WAVI - Song Meter Micro
wavi_calltypes_sm <- c("Song", "Song", "Call", "Song", "Song")
names(wavi_calltypes_sm) <- 1:5
wavi_sm$calltype <- as.character(wavi_calltypes_sm[wavi_sm$cluster])



# adds a datetime string (YYYYMMDD_HHMMSS) to dataframe 
addDateTimeLoc <- function(metadata)
{
  metadata$datetimeloc <- paste(metadata$site,
                             "_",
                             sprintf("%04d", metadata$year),
                             sprintf("%02d", metadata$month),
                             sprintf("%02d", metadata$day),
                             "_",
                             sprintf("%02d", metadata$hour),
                             sprintf("%02d", metadata$minute),
                             sprintf("%02d", metadata$second), 
                             sep="")
  return(metadata)
}

# Code to reconstruct filepaths to audiofiles based on metadata in dataframe
#   this time with output files for full 1 minute clips
getAudioFilepathFinal <- function(metadata, audio_dir, common_name)
{
  # Construct date and time for filename from metadata
  filename_bases <- paste(sprintf("%04d", metadata$year),
                          sprintf("%02d", metadata$month),
                          sprintf("%02d", metadata$day),
                          "_",
                          sprintf("%02d", metadata$hour),
                          sprintf("%02d", metadata$minute),
                          sprintf("%02d", metadata$second), 
                          sep="")
  # Double-check that year is correct in audio path
  audio_dir <- str_replace_all(audio_dir, "2022", as.character(metadata$year))
  audio_dir <- str_replace_all(audio_dir, "2023", as.character(metadata$year))
  # Append datetime to directory structure to reconstruct audio filepath
  original_filepaths <- paste(audio_dir, "/", metadata$site, "/", filename_bases, ".wav", sep="")
  original_detection_filepaths <- paste(audio_dir, "/../detections/", metadata$site, "/", filename_bases, ".BirdNET.results.csv", sep="")
  # Create new filename to save test 3-sec audio snippet to 
  new_filepaths <- paste("G:/Bioacoustics/Goleta_2023/training_data/", common_name, "/audio/", metadata$unit_type, "_cluster_", metadata$cluster, "_", metadata$site, "_", filename_bases, "_", metadata$start_s, ".wav", sep="")
  new_detection_filepaths <- paste("G:/Bioacoustics/Goleta_2023/training_data/", common_name, "/detections/", metadata$unit_type, "_cluster_", metadata$cluster, "_", metadata$site, "_", filename_bases, "_", metadata$start_s, ".csv", sep="")
  
  return(data.frame(new_filepath = new_filepaths,
                    filename_base = filename_bases,
                    original_filepath = original_filepaths,
                    new_detection_filepath = new_detection_filepaths,
                    original_detection_filepath = original_detection_filepaths,
                    confidence = metadata$confidence,
                    cluster = metadata$cluster,
                    start_s = metadata$start_s,
                    unit_type = metadata$unit_type,
                    common_name = common_name))
}

# For a given bird species and a given unique datetime
#   Aggregate all calls of a given cluster 
#   Count how many occur of each type 
extractTrainingAudio <- function(embeddings_dataset, species, total_samples=100)
{
  set.seed(1)
  output_datetimelocs <<- list()
  
  ## Double-check whether output directory exists, and if not, create it:
  if(!(dir.exists(paste(audio_directory, "/..", "/training_data", "/", species, "/audio/", sep=""))))
    dir.create(paste(audio_directory, "/..", "/training_data", "/", species, "/audio/", sep=""), 
               showWarnings = FALSE, recursive = TRUE)
  if(!(dir.exists(paste(audio_directory, "/..", "/training_data", "/", species, "/detections/", sep=""))))
    dir.create(paste(audio_directory, "/..", "/training_data", "/", species, "/detections/", sep=""), 
               showWarnings = FALSE, recursive = TRUE)
  
  # Add unique audio file identifiers (sitename, date, and timestamp)
  embeddings_dataset <- addDateTimeLoc(embeddings_dataset)
  # Get list of all unique audiofiles which contained the target species
  all_datetimelocs <- unique(embeddings_dataset$datetimeloc)
  # Length of the above is the max number of audio files which can be extracted
  max_files <- length(all_datetimelocs)
  
  # Get number of unique vocal cluster types
  num_clusters <- length(unique(embeddings_dataset$cluster))
  
  # Iterate over each cluster type
  for(cluster_ind in 1:length(unique(embeddings_dataset$cluster)))
  {
    print(paste("  Beginning to work on cluster ", cluster_ind, sep=""))
    cluster_count <<- 0
    # Filter to target cluster
    cluster_data <- embeddings_dataset %>% filter(cluster == cluster_ind)
    # Code to gather example audio files for a given confidence range and cluster type
    iterativelyGatherExamples <- function(dataset, quantile_low, quantile_high, num_subset_samples)
    {
      current_samples <- 0
      while((current_samples < num_subset_samples))
      {
        target_dataset <- cluster_data %>% filter(!(datetimeloc %in% output_datetimelocs),
                                               confidence <= quantile_high,
                                               confidence > quantile_low)
        # Escape if no good data remains in this category 
        #   (this might happen for example if there are no high-confidence incidences in a given call type)
        if(nrow(target_dataset) == 0)
          break;
        # Choose one random bird detection from appropriate data
        target_datapoint <- target_dataset[sample(1:nrow(target_dataset),1),]
        # Get filepath to source audio and filepath to copy audio to
        target_file_info <- getAudioFilepathFinal(target_datapoint, audio_directory, species)
        file.copy(target_file_info$original_filepath, target_file_info$new_filepath)
        file.copy(target_file_info$original_detection_filepath, target_file_info$new_detection_filepath)
        # Increment counters
        current_samples <- current_samples+1
        output_datetimelocs <<- c(output_datetimelocs, target_datapoint$datetimeloc)
        cluster_count <<- cluster_count + 1
      }
    }
    
    # Start from the highest confidence quantiles (0.9 to 1.0)
    # Iteratively add clips from each bin, trying to get total_samples/9 for each range
    # For each subsequent lower range, add extra requested samples to accomodate those which were not retrieved for higher ranges
    # This is because for scarce species, there may be too few high-confidence detections to fill a request
    # This will ensure that the overall number of samples is total_samples, while also ensuring that
    #   the spread is as even as possible across confidence bins
    iterativelyGatherExamples(cluster_data, 0.9, 1.0, total_samples/num_clusters/9)
    print(paste("    Number of songs after confidence range 0.9 to 1.0: ", cluster_count, sep=""))
    iterativelyGatherExamples(cluster_data, 0.8, 0.9, total_samples/num_clusters/9 + ((total_samples/num_clusters/9)-cluster_count)/8)
    print(paste("    Number of songs after confidence range 0.8 to 0.9: ", cluster_count, sep=""))
    iterativelyGatherExamples(cluster_data, 0.7, 0.8, total_samples/num_clusters/9 + ((total_samples/num_clusters*2/9)-cluster_count)/7)
    print(paste("    Number of songs after confidence range 0.7 to 0.8: ", cluster_count, sep=""))
    iterativelyGatherExamples(cluster_data, 0.6, 0.7, total_samples/num_clusters/9 + ((total_samples/num_clusters*3/9)-cluster_count)/6)
    print(paste("    Number of songs after confidence range 0.6 to 0.7: ", cluster_count, sep=""))
    iterativelyGatherExamples(cluster_data, 0.5, 0.6, total_samples/num_clusters/9 + ((total_samples/num_clusters*4/9)-cluster_count)/5)
    print(paste("    Number of songs after confidence range 0.5 to 0.6: ", cluster_count, sep=""))
    iterativelyGatherExamples(cluster_data, 0.4, 0.5, total_samples/num_clusters/9 + ((total_samples/num_clusters*5/9)-cluster_count)/4)
    print(paste("    Number of songs after confidence range 0.4 to 0.5: ", cluster_count, sep=""))
    iterativelyGatherExamples(cluster_data, 0.3, 0.4, total_samples/num_clusters/9 + ((total_samples/num_clusters*6/9)-cluster_count)/3)
    print(paste("    Number of songs after confidence range 0.3 to 0.4: ", cluster_count, sep=""))
    iterativelyGatherExamples(cluster_data, 0.2, 0.3, total_samples/num_clusters/9 + ((total_samples/num_clusters*7/9)-cluster_count)/2)
    print(paste("    Number of songs after confidence range 0.2 to 0.3: ", cluster_count, sep=""))
    iterativelyGatherExamples(cluster_data, 0.0, 0.2, total_samples/num_clusters - cluster_count)
    print(paste("    Number of songs after confidence range 0.0 to 0.2: ", cluster_count, sep=""))
  }
}

# Generate training subsets for each riparian species of interest
extractTrainingAudio(yewa_am, "yellow_warbler", 500)
extractTrainingAudio(yewa_sm, "yellow_warbler", 500)
extractTrainingAudio(wiwa_am, "wilsons_warbler", 500)
extractTrainingAudio(wiwa_sm, "wilsons_warbler", 500)
extractTrainingAudio(wavi_am, "warbling_vireo", 500)
extractTrainingAudio(wavi_sm, "warbling_vireo", 500)
extractTrainingAudio(cbch_am, "chestnut_backed_chickadee", 500)
extractTrainingAudio(cbch_sm, "chestnut_backed_chickadee", 500)
extractTrainingAudio(pufi_am, "purple_finch", 500)
extractTrainingAudio(pufi_sm, "purple_finch", 500)
extractTrainingAudio(sosp_am, "song_sparrow", 500)
extractTrainingAudio(sosp_sm, "song_sparrow", 500)
extractTrainingAudio(acwo_am, "acorn_woodpecker", 500)
extractTrainingAudio(acwo_sm, "acorn_woodpecker", 500)
extractTrainingAudio(ocwa_am, "orange_crowned_warbler", 500)
extractTrainingAudio(ocwa_sm, "orange_crowned_warbler", 500)
extractTrainingAudio(spto_am, "spotted_towhee", 500)
extractTrainingAudio(spto_sm, "spotted_towhee", 500)



# Look up embedding for a given detection
lookUpEmbedding <- function(metadata, return_embedding=FALSE)
{
  # fill this helper function out later
}

# Check over all detections in sample
#   NOTE this includes ALL detections, not just the target species
#   Currently, this assigns all detections the cluster of the detection used to select this audio segment! 
all_yewa_training_data <- bind_rows(lapply(list.files("G:/Bioacoustics/Goleta_2023/training_data/yellow_warbler/detections/",
                                                      full.names=TRUE), 
                                           function(filename) {
                                             new_data <- read_csv(filename, col_types=cols()) %>%
                                               janitor::clean_names()
                                             new_data$filename <- filename
                                             new_data$site <- substr(basename(filename), 14, 15)
                                             new_data$unit_type <- substr(basename(filename), 1, 2)
                                             return(new_data)
                                           })) %>%
  lookUpEmbedding(return_embedding=FALSE)
# Visualize distribution of confidence levels in training set
ggplot(all_yewa_training_data %>% filter(common_name == "Yellow Warbler")) + 
  geom_histogram(aes(x=confidence)) + 
  facet_wrap(~cluster)
# Visualize number of detections (any species) / recording in training set 
ggplot(all_yewa_training_data %>% group_by(filename, cluster) %>% tally()) + 
  geom_histogram(aes(x=n)) + 
  facet_wrap(~cluster)
# Visualize number of detections (target species) / recording in training set 
ggplot(all_yewa_training_data %>% filter(common_name =="Yellow Warbler") %>% group_by(filename, cluster) %>% tally()) + 
  geom_histogram(aes(x=n)) + 
  facet_wrap(~cluster)

# Extract a subset of target dataset into a folder for someone to look at and manually validate
existing_audio_dir <- "G:/Bioacoustics/Goleta_2023/training_data/initial_test_set/"
output_audio_dir <- "G:/Bioacoustics/Goleta_2023/training_data/zach_data_assignment"
num_samples <- 100 # for each species 
species_list <- list.files(existing_audio_dir, full.names=TRUE)

# Move [num_samples] random paired audio/detections files from each species into assignment folder
#   Removes these files from the source directory so they're only assigned once
#   Note it's theoretically possible (although unlikely) that some files could match across species
#   In that case only one copy will end up getting transfered (which makes sense, it only needs to be listened to once)
lapply(species_list,
       function(target_directory){
         # Double-check that target is a species; otherwise continue
         if(dir.exists(target_directory))
         {
           # Move audio data
           audio_files <- list.files(paste(target_directory, "/audio/", sep=""))
           random_indices <- sample(1:length(audio_files), num_samples)
           original_filenames <- paste(target_directory, "/audio/", audio_files[random_indices], sep="")
           new_filenames <- paste(output_audio_dir, "/audio/", audio_files[random_indices], sep="")
           file_successes <- file.rename(original_filenames, new_filenames)
           # Move detections data
           detection_files <- list.files(paste(target_directory, "/detections/", sep=""))
           original_detection_filenames <- paste(target_directory, "/detections/", detection_files[random_indices], sep="")
           new_detection_filenames <- paste(output_audio_dir, "/detections/", detection_files[random_indices], sep="")
           detections_successes <- file.rename(original_detection_filenames, new_detection_filenames)
           # Exit, indicating success
           return(1)
         }
         return(-1)
       })

# Add extra fields to detections files for volunteers to fill out, and remove extra fields
# Retains:
#   start_s           (start time of period in which detection occurred)
#   end_s             (end time of period in which detection occurred)
#   common_name       (common name of bird species)
#   correct           (user input - is species correct?)
#   true_common_name  (user input - if NOT correct, what is actual species?)
#   sound_type        (user input - either Call or Song)
addValidationFields <- function(filename)
{
  # Read input data
  detections_data <- read_csv(filename, col_types=cols()) %>%
    # standardize column heading format
    janitor::clean_names() %>% 
    # remove extra columns
    select(c("start_s", "end_s", "common_name")) 
  # Add user input columns
  detections_data$correct <- ""
  detections_data$true_common_name <- ""
  detections_data$sound_type <- ""
  # Write output 
  write_csv(detections_data, filename)
}

# Apply the function 
training_detections_set <- lapply(list.files(paste(output_audio_dir, "/detections", sep=""), full.names=TRUE),
                                  addValidationFields)


