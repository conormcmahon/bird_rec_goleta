
# Dig in fractal manner through tree of directories to a target depth (in directories)
# Then split all .WAV files at terminal level to at most 1 minute duration 

# DIRECTORY - base level directory for searches
# LEVELS_TO_TRAVERSE - number of directories to delve down through to level with audio files
# CURRENT_LEVEL - start off at 1, automatically increments within recursive function below
traverseDirectories <- function(directory, levels_to_traverse, current_level=1)
{
  # Get files and directories at this level
  current_files <- list.files(directory)
  current_dirs <- list.dirs(directory, recursive=FALSE)
  print("")
  print(paste(rep("  ", current_level), "Current working direcotry:", sep=""))
  print(directory)
  print(paste(rep("  ", current_level), "  Current level is ", current_level, " of ", levels_to_traverse, sep=""))
  print(paste(rep("  ", current_level), "  Number of files: ", length(current_files), sep=""))
  print(paste(rep("  ", current_level), "  Number of directories: ", length(current_dirs), sep=""))
  
  # If we've reached terminal level, apply splitting function
  if(levels_to_traverse == current_level)
    return(lapply(paste(directory, "/", current_files, sep=""), splitAudioFile, length_seconds=60))
  # Otherwise, dig to another level...
  return(lapply(current_dirs, traverseDirectories, levels_to_traverse=levels_to_traverse, current_level=current_level+1))
}

results_2022 <- traverseDirectories("G:/Bioacoustics/Goleta_2022/audio/", 2, 1)
results_2023 <- traverseDirectories("G:/Bioacoustics/Goleta_2023/audio/", 2, 1)
