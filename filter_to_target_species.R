
# Filters a list of detections to include only birds from a given target list 

library(tidyverse)
library(janitor)
library(sf)

# Superdirectory containing directories for each checklist
data_directory <- "G:/Bioacoustics/SCR_2023/detections/"
output_directory <- "G:/Bioacoustics/SCR_2023/detections_riparian/"

# List of all target classes
target_species <- c("Sayornis nigricans",
                    "Melospiza melodia",
                    "Pheucticus melanocephalus",
                    "Pipilo maculatus",
                    "Contopus sordidulus",
                    "Molothrus ater",
                    "Icterus bullockii",
                    "Empidonax difficilis",
                    "Geothlypis trichas",
                    "Vireo gilvus",
                    "Cardellina pusilla",
                    "Tachycineta bicolor",
                    "Setophaga petechia",
                    "Haemorhous purpureus",
                    "Vireo bellii",
                    "Coccyzus americanus",
                    "Catharus ustulatus",
                    "Passerina caerulea",
                    "Empidonax traillii",
                    "Icteria virens",
                    "Archilochus alexandri",
                    "Leuconotopicus villosus",
                    "Poecile rufescens",
                    "Catherpes mexicanus")

# Ingest a single recording from BirdNet .csv tab-delimited report
outputFilteredDetections <- function(filename, directory, output_directory, checklist_name)
{
  data <- read_csv(paste(directory,filename,sep=""),
                   show_col_types = FALSE) 
  data <- data[which(unlist(lapply(as.vector(data[,3]), function(name){ name %in% target_species}))),]
  data$filename <- rep(paste(directory,filename,sep="/"), nrow(data))
  data$site <- rep(checklist_name, nrow(data))
  write_csv(data, paste(output_directory,filename,sep=""))
  return(data)
}
# Process all recording files associated with a single checklist directory
filterDetections <- function(checklist)
{
  checklist_directory <- paste(data_directory,checklist,"/",sep="")
  files <- list.files(path=checklist_directory)
  final_directory <- paste(output_directory,checklist,"/",sep="")
  csv_files <- files[grepl(".*\\.csv", files)] 
  dir.create(file.path(output_directory, checklist), showWarnings=FALSE)
  
  dataframe_list <- lapply(csv_files, outputFilteredDetections, directory=checklist_directory, output_directory=final_directory, checklist_name=checklist)
  # Drop any empty dataframes (no birds detected)
  dataframe_list <- dataframe_list[lapply(dataframe_list, nrow) > 0]
  # Convert list of dataframes to one dataframe
  dataframe_out <- bind_rows(dataframe_list)
  return(dataframe_out)
}
# Get all checklist IDs within data_directory
checklists <- list.files(path=data_directory)
# Drop any .csv readme or config files
checklists <- checklists[(1:length(checklists))*(1-grepl(".*\\.csv", checklists))]
# Collect all data from files
all_detections <- lapply(checklists, filterDetections)
# Remove any empty checklists
all_detections <- bind_rows(all_detections)

