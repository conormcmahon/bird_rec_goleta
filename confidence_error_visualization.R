
# Load labelled data from manual listening validation
# Perform some simple sanity checks and summary statistics to understand how confidence and correctness are related by species

library(tidyverse)

# Load detections data, including confidence value
all_detections <- rbind(read_csv("G:/Bioacoustics/Goleta_2022/detections/all_detections_compiled.csv"),
                        read_csv("G:/Bioacoustics/Goleta_2023/detections/all_detections_compiled_2023.csv"))
# Load labelled data, which do not include confidences
labelled_data <- read_csv("G:/Bioacoustics/Goleta_2023/labelled_data/labelled_call_data.csv")

# Add site and time metadata to labelled dataset, extracted from filenames
lookupInfo <- function(filename)
{
  # Break filename by _ character
  filename_chunks <- strsplit(as.character(filename), "_")[[1]]
  # Get date and time from filename
  year_c <- as.numeric(substr(filename_chunks[6],1,4))
  month_c <- as.numeric(substr(filename_chunks[6],5,6))
  day_c <- as.numeric(substr(filename_chunks[6],7,8))
  hour_c <- as.numeric(substr(filename_chunks[7],1,2))
  minute_c <- as.numeric(substr(filename_chunks[7],3,4))
  second_c <- as.numeric(substr(filename_chunks[7],5,6))
  # Combine all data into new dataframe
  metadata <- data.frame(site = paste(filename_chunks[4], "_", filename_chunks[5], sep=""),
                         year = year_c,
                         month = month_c,
                         day = day_c,
                         hour = hour_c,
                         minute = minute_c,
                         second = second_c)
  return(metadata)
}
# Run the above function to get metadata for each labelled observation, merging back with original dataframe
labelled_data <- cbind(labelled_data, 
                       bind_rows(lapply(1:nrow(labelled_data), 
                                        function(ind){ 
                                          return(lookupInfo(labelled_data[ind,"filename"]))
                                          })))

# Merge the above datasets to get labels and confidence together
labelled_data <- merge((labelled_data %>% 
                          dplyr::select("common_name", "site",
                                        "year", "month", "day", 
                                        "hour", "minute", "second", 
                                        "start_s", "correct", 
                                        "true_common_name", "sound_type", "primary_type")), 
                       all_detections)


# Create a table showing each species and number of correct or incorrect detections
labelled_data %>%
  filter(primary_type %in% c("Call", "Song"))  %>%
  group_by(common_name,
           (correct > 0)) %>%
  tally() %>%
  pivot_wider(names_from='(correct > 0)', values_from=n) %>% 
  View("Error by Species")

# Create a table showing the incorrect detections which had the highest-confidence 
labelled_data %>% 
  dplyr::select(-c("sound_type", "time_frac")) %>%
  filter(correct != 1) %>% 
  arrange(-confidence) %>% 
  View("Worst Errors")

# Now, for focal species, visualize breakdown of 'correctness' of individual calls vs. BirdNET confidence:
#   For a given species, create plot and logistic regression of confidence vs. correctness
visualizeConfidence <- function(target_species)
{
  # Filter dataset to species
  target_data <- labelled_data %>% 
    filter(common_name == target_species,
           !is.na(correct))
  # Perform logistic regression
  sp_logit <- glm(data=target_data,
                    (correct == 1) ~ confidence)
  print(sp_logit)
  print(summary(sp_logit))
  # Generate predictions from logit model to visualize range of expected likelihood-confidence relationships
  #   Force model to compute standard errors for each value, as well
  prediction <- predict(sp_logit, 
                          newdata = data.frame(confidence = (1:99)/100),
                          se=TRUE)
  prediction_df <- data.frame(confidence = (1:99)/100, 
                              fit = prediction$fit,
                              se = prediction$se.fit)
  # Create plot 
  sp_plot <- ggplot() + 
    geom_point(data = target_data,
               aes(x=confidence, y=(correct == 1)*1)) + 
    geom_line(data = prediction_df,                        # Logit prediction
              aes(x=confidence, y=fit)) + 
    geom_line(data = prediction_df,                        # Upper 95th percentile
              aes(x=confidence, y=fit+1.96*se), col="red") + 
    geom_line(data = prediction_df,                        # Lower 95th percentile
              aes(x=confidence, y=fit-1.96*se), col="red") + 
    scale_y_continuous(limits=c(0,1)) + 
    ggtitle(paste("Confidence Regression for ", target_species, sep=""))
  print(sp_plot)
  
  return(list(sp_logit,
              prediction_df,
              sp_plot))
}

# Get list of all errors for a given species
getSpeciesErrors <- function(target_species)
{
  labelled_data %>% 
    dplyr::select(-c("sound_type", "time_frac", "audio_chunk")) %>%
    filter(correct != 1) %>% 
    arrange(-confidence) %>% 
    filter(common_name == target_species)
}
# Get summary of correct vs. incorrect detections for a given species
getSpeciesErrorSummary <- function(target_species)
{
  labelled_data %>%
    filter(primary_type %in% c("Call", "Song"),
           common_name == target_species)  %>%
    group_by((correct > 0)) %>%
    tally() %>%
    pivot_wider(names_from='(correct > 0)', values_from=n)
}


# ************** Summary by Species **************

# **** Orange-crowned Warbler ****
target_species <- "Orange-crowned Warbler"
getSpeciesErrorSummary(target_species)
#   15 errors / 229 detections --> 6.55%
getSpeciesErrors(target_species)
#   worst error at 0.488 confidence
ocwa_results <- visualizeConfidence(target_species)


# Acorn Woodpecker
target_species <- "Acorn Woodpecker"
getSpeciesErrorSummary(target_species)
#   0 errors / 139 detections --> 0.00%
getSpeciesErrors(target_species)
#   no errors so far
acwo_results <- visualizeConfidence(target_species)


# Spotted Towhee
target_species <- "Spotted Towhee"
getSpeciesErrorSummary(target_species)
#   0 errors / 263 detections --> 0.00%
getSpeciesErrors(target_species)
#   no errors so far
spto_results <- visualizeConfidence(target_species)


# Wilson's Warbler
target_species <- "Wilson's Warbler"
getSpeciesErrorSummary(target_species)
#   11 errors / 119 detections --> 9.24%
getSpeciesErrors(target_species)
#   worst error at 0.689 confidence
wiwa_results <- visualizeConfidence(target_species)


# Warbling Vireo
target_species <- "Warbling Vireo"
getSpeciesErrorSummary(target_species)
#   6 errors / 80 detections --> 7.5%
getSpeciesErrors(target_species)
#   worst error at 0.377 confidence
wavi_results <- visualizeConfidence(target_species)


# Yellow Warbler
target_species <- "Yellow Warbler"
getSpeciesErrorSummary(target_species)
#   0 errors / 80 detections --> 0.00%
getSpeciesErrors(target_species)
#   no errors so far
yewa_results <- visualizeConfidence(target_species)

