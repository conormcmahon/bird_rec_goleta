
library(tidyverse)
library(here)
library(lubridate)
library(janitor)

yewa <- read_csv(here::here("bbs_data","yewa.csv"))
wiwa <- read_csv(here::here("bbs_data","wiwa.csv"))
wavi <- read_csv(here::here("bbs_data","wavi.csv"))
casj <- read_csv(here::here("bbs_data","casj.csv"))
calt <- read_csv(here::here("bbs_data","calt.csv"))
caqu <- read_csv(here::here("bbs_data","caqu.csv"))
wren <- read_csv(here::here("bbs_data","wren.csv"))
sosp <- read_csv(here::here("bbs_data","sosp.csv"))

# Add day of year to dataframe
addDOY <- function(data_in)
{
  dates_list <- str_split(data_in$observation_date, " ")
  getIthElement <- function(data_vector, index)
  {
    return(data_vector[[index]])
  }
  data_in$day <- unlist(lapply(dates_list, getIthElement, index=1))
  data_in$month_alph <- unlist(lapply(dates_list, getIthElement, index=2))
  months <- 1:12
  names(months) <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  data_in$month <- months[data_in$month_alph]
  data_in$year <- unlist(lapply(dates_list, getIthElement, index=3))
  data_in$date <- as.Date(paste(data_in$year, data_in$month, data_in$day, sep="-"))
  data_in$doy <- lubridate::yday(data_in$date)
  
  return(data_in)
}

all_species <- bind_rows(yewa, wiwa, wavi, casj, calt, caqu, wren, sosp) %>%
  janitor::clean_names()
all_species <- addDOY(all_species)

ggplot() + 
  geom_histogram(data=all_species %>% filter(code %in% c("YEWA","WIWA","WAVI")),
                 aes(x=doy)) + 
  facet_wrap(~common_name) + 
  labs(x="Day of Year",
       y="Breeding Record Counts",
       title="Seasonality of Riparian Breeding Bird Records",
       )
ggplot() + 
  geom_histogram(data=all_species %>% filter(code %in% c("SOSP","CASJ","CALT","WREN","CAQU")),
                 aes(x=doy)) + 
  facet_wrap(~common_name) + 
  labs(x="Day of Year",
       y="Breeding Record Counts",
       title="Seasonality of Riparian Breeding Bird Records",
  )
  

