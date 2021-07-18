
library(tidyverse)
library(janitor)

# Superdirectory containing directories for each checklist
data_directory <- "D:/bird_rec/data/2021_03_san_jose_creek/"

# Ingest a single recording from BirdNet .txt tab-delimited report
ingestRecording <- function(filename, directory)
{
  data <- read_tsv(paste(directory,filename,sep="/"),col_types=cols()) %>% janitor::clean_names()
  return(data)
}
# Process all recording files associated with a single checklist directory
getBirdData <- function(checklist)
{
  checklist_directory <- paste(data_directory,checklist,sep="/")
  files <- list.files(path=checklist_directory)
  text_files <- files[grepl(".*\\.txt", files)] 
  
  dataframe_list <- lapply(text_files, ingestRecording, directory=checklist_directory)
  # Drop any empty dataframes (no birds detected)
  dataframe_list <- dataframe_list[lapply(dataframe_list, nrow) > 0]
  # Convert list of dataframes to one dataframe
  bind_rows(dataframe_list)
}
# Get all checklist IDs within data_directory
checklists <- list.files(path=data_directory)
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
              best_rank=min(rank),
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
            presence_rate=n()/length(checklists)) %>% 
  arrange(-total)
view(species_presence)
