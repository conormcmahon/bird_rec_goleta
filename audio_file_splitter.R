
# Much of this code based on the helpful StackExchange thread by @Jota, here:
#    https://stackoverflow.com/questions/20713513/efficiently-split-a-large-audio-file-in-r

library(seewave)
library(tuneR)
library(lubridate)
library(stringr)

splitAudioFile <- function(filepath, length_seconds)
{
  # Check that target file ends in .WAV extension
  if(!(substr(filepath, nchar(filepath)-3, nchar(filepath)) %in% c(".wav", ".Wav", ".WAV")))
    return(-1)
  
  audio <- readWave(filepath)
  
  # the frequency of your audio file
  freq <- audio@samp.rate
  # the length and duration of your audio file
  totlen <- length(audio)
  totsec <- totlen/freq 
  
  # If target file is less than target length, return without changing it
  #   UNLESS we need to clean up the name 
  if(totsec <= length_seconds) 
  {
    # Check whether there's any messy prefix on the filename before the datetime
    filename_base <- substr(filepath, 1, nchar(filepath)-19)
    # Slice off any information which comes BEFORE the datetime
    folder_cutoff_character_ind <- str_locate_all(pattern="/", filename_base)
    folder_cutoff_character_ind <- folder_cutoff_character_ind[[1]][[nrow(folder_cutoff_character_ind[[1]]), 1]]
    if(nchar(filename_base) == nchar(substr(filename_base, 1, folder_cutoff_character_ind)))
      return(0)
    
    # Save the file with new clipped name
    filename_out <- paste(substr(filename_base, 1, folder_cutoff_character_ind),
                          substr(filepath, nchar(filepath)-18, nchar(filepath)), sep="")
    writeWave(audio, 
              filename=filename_out)
    return(0)
  }
  
  # defining the break points
  breaks <- unique(c(seq(0, totsec, length_seconds), totsec))
  index <- 1:(length(breaks)-1)
  
  # the split
  leftmat<-matrix(audio@left, ncol=(length(breaks)-2), nrow=length_seconds*freq) 
  # the warnings are nothing to worry about here... 
  
  # convert to list of Wave objects.
  subsamps0409_180629 <- lapply(1:ncol(leftmat), function(x)Wave(left=leftmat[,x], 
                                                                 samp.rate=audio@samp.rate,
                                                                 bit=audio@bit)) 
  
  
  # get the last part of the audio file.  the part that is < length_seconds
  lastbit <- audio@left[(breaks[length(breaks)-1]*freq):length(audio)]
  
  # convert and add the last bit to the list of Wave objects
  subsamps0409_180629[[length(subsamps0409_180629)+1]] <- 
    Wave(left=lastbit, samp.rate=audio@samp.rate, bit=audio@bit)
  
  # finally, save the Wave objects
  #setwd("C:/Users/Whatever/Wave_object_folder")
  
  # I had some memory management issues on my computer when doing this
  # process with large (~ 130-150 MB) audio files so I used rm() and gc(),
  # which seemed to resolve the problems I had with allocating memory.
  rm("audio","freq","index","lastbit","leftmat","totlen","totsec")
  
  gc()
  
  # Set up new filenames using new timestamps for each subset
  # Base filename from input
  filename_base <- substr(filepath, 1, nchar(filepath)-19)
  # Slice off any information which comes BEFORE the datetime
  folder_cutoff_character_ind <- str_locate_all(pattern="/", filename_base)
  folder_cutoff_character_ind <- folder_cutoff_character_ind[[1]][[nrow(folder_cutoff_character_ind[[1]]), 1]]
  filename_base <- substr(filename_base, 1, folder_cutoff_character_ind)
  # Timestamp and date of start file
  timestamp_start <- substr(filepath, nchar(filepath)-9, nchar(filepath)-4)
  date_start <- substr(filepath, nchar(filepath)-18, nchar(filepath)-11)
  datetime_start = lubridate::ymd_hms(paste(substr(date_start, 1, 4), "-",
                                            substr(date_start, 5, 6), "-",
                                            substr(date_start, 7, 8), " ",
                                            substr(timestamp_start, 1, 2), ":",
                                            substr(timestamp_start, 3, 4), ":",
                                            substr(timestamp_start, 5, 6), ":"))
  # New vector of filenames, including updated timestamps
  #   sprintf, when used in this format, forces output to use leading zeroes where appropriate
  #   e.g. Feb. 5, 2023 will be 20230205, rather than 202325
  file_datetimes <- datetime_start + as.integer(round(breaks))
  filenames <- paste(filename_base,
                     sprintf("%02d", year(file_datetimes)),
                     sprintf("%02d", month(file_datetimes)),
                     sprintf("%02d", day(file_datetimes)), "_",
                     sprintf("%02d", hour(file_datetimes)),
                     sprintf("%02d", minute(file_datetimes)),
                     sprintf("%02d", second(file_datetimes)),
                     ".WAV",sep="")
  
  # Save the files
  sapply(1:length(subsamps0409_180629),
         function(x)writeWave(subsamps0409_180629[[x]], 
                              filename=filenames[x]))
  
  return(length(breaks))
  
  rm(subsamps0409_180629)
  gc()
}