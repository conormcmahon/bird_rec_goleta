
library(tidyverse)
library(here)
library(janitor)

# User-defined parameters for selecting likely true IDs
confidence_threshold <- 0.85
detections_threshold <- 250

# Load list of bird names
species_metadata <- read_csv("D:/birdrec/scripts/birdnet_species_list_processing/all_species_metadata.csv")
bird_names <- (species_metadata %>% filter(bird))$common_name

# Load a list which includes the species from the Riparian Broad species list
birdnet_riparian_broad_list <- read_csv("D:/birdrec/protocol/birdnet_riparian_broad.csv")

# Load another list which includes commonly-reported birds from eBird in the county which are missing above
missing_from_birdnet_list <- read_csv("D:/birdrec/protocol/missing_from_birdnet_riparian_broad.csv") %>%
  filter(num_reports > 2500)

# Combine the above into an 'allowed species' list
allowed_species <- unique(c(birdnet_riparian_broad_list$common_name,
                            missing_from_birdnet_list$common_name,
                            "Mountain Quail",
                            "Greater Roadrunner",
                            "Cassin's Vireo",
                            "Olive-sided Flycatcher"))

# Load species guild designations
species_guilds <- read_csv("D:/birdrec/protocol/species_guilds.csv")
species_guilds <- species_guilds %>% 
  replace(is.na(.), 0) %>%
  janitor::clean_names()

# Load effort data
effort <- rbind(read_csv("G:/Bioacoustics/Goleta_2022/survey_time.csv"),
                read_csv("G:/Bioacoustics/Goleta_2023/survey_time.csv"),
                read_csv("G:/Bioacoustics/Goleta_2024/survey_time.csv"))
effort_daily_summary <- effort %>%
  group_by(site) %>% 
  summarize(hours_surveyed = n()/60)

# Load detections
all_detections <- rbind(read_csv("G:/Bioacoustics/Goleta_2022/detections_compiled/all_detections.csv"),
                        read_csv("G:/Bioacoustics/Goleta_2023/detections_compiled/all_detections.csv"),
                        read_csv("G:/Bioacoustics/Goleta_2024/detections_compiled/all_detections.csv")) %>%
  filter(common_name %in% allowed_species)

# Select just Ellwood, Atascadero, Rattlesnake, and Mission Creek data
ew_ac_detections <- all_detections %>% 
  filter(substr(site, 1,2) %in% c("EW","AC","RC","MC"))
# Add dominant vegetation type to each 
veg_types <- c("willow", "eucalyptus", "oak-laurel", "oak-laurel")
names(veg_types) <- c("AC", "EW", "RC", "MC")
ew_ac_detections$vegetation <- veg_types[substr(ew_ac_detections$site,1,2)]

# Get species feeding guilds
ggplot(ew_ac_detections) + 
  geom_histogram(aes(x=yday)) + 
  facet_wrap(~(substr(site,1,2)=="EW"))


# Get daily summary
hourly_detections <- ew_ac_detections %>% 
  group_by(site, year, yday, hour, common_name) %>%
  summarize(count = n(), 
            conf_100 = max(confidence),
            conf_90 = quantile(confidence, 0.9),
            conf_90 = quantile(confidence, 0.75),
            conf_90 = quantile(confidence, 0.5),
            conf_90 = quantile(confidence, 0.25),
            conf_mean = mean(confidence))
hourly_effort <- effort %>% 
  mutate(yday = doy) %>%
  select(-doy) %>%
  group_by(site, year, yday, hour) %>%
  summarize(recording_count = n(),
            time_recorded = sum(length))

# Get hourly density for a target species
getHourlyDensity <- function(target_species)
{
  species_density <- merge(hourly_detections %>% filter(common_name == target_species),
                           hourly_effort, 
                           by=c("site","year","yday","hour"),
                           all=TRUE) %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    mutate(call_density = count / recording_count) %>%
    mutate(common_name = target_species)
  species_density$vegetation <- veg_types[substr(species_density$site,1,2)]
  return(species_density)
}

yewa_density <- getHourlyDensity("Yellow Warbler")
ocwa_density <- getHourlyDensity("Orange-crowned Warbler")
sosp_density <- getHourlyDensity("Song Sparrow")
bhgr_density <- getHourlyDensity("Black-headed Grosbeak")
psfl_density <- getHourlyDensity("Pacific-slope Flycatcher")
casj_density <- getHourlyDensity("California Scrub-Jay")
hofi_density <- getHourlyDensity("House Finch")
lego_density <- getHourlyDensity("Lesser Goldfinch")
spto_density <- getHourlyDensity("Spotted Towhee")
hoor_density <- getHourlyDensity("Hooded Oriole")
anhu_density <- getHourlyDensity("Anna's Hummingbird")
alhu_density <- getHourlyDensity("Allen's Hummingbird")


# Orange-crowned Warbler
summary(lm(data = ocwa_density %>% 
             filter(year==2024), 
           call_density ~ sin((hour+3)*pi/24) * yday + vegetation))
ggplot(ocwa_density %>% filter(year==2024)) + 
  geom_point(aes(x=hour, y=call_density), alpha=0.1) + 
  geom_smooth(aes(x=hour, y=call_density), span = 0.1) + 
  facet_wrap(~vegetation)

# Pacific-slope Flycatcher
summary(lm(data = psfl_density %>% 
             filter(year==2024), 
           call_density ~ sin((hour+3)*pi/24) * yday + vegetation))
ggplot(psfl_density %>% filter(year==2024)) + 
  geom_point(aes(x=hour, y=call_density), alpha=0.1) + 
  geom_smooth(aes(x=hour, y=call_density), span = 0.1) + 
  facet_wrap(~vegetation)

# Song Sparrow
summary(lm(data = sosp_density %>% 
             filter(year==2024), 
           call_density ~ sin((hour+3)*pi/24) * yday + vegetation))
ggplot(sosp_density %>% filter(year==2024)) + 
  geom_point(aes(x=hour, y=call_density), alpha=0.1) + 
  geom_smooth(aes(x=hour, y=call_density), span = 0.1) + 
  facet_wrap(~vegetation)

# House Finch
summary(lm(data = hofi_density %>% 
             filter(year==2024), 
           call_density ~ sin((hour+3)*pi/24) * yday + vegetation))
ggplot(hofi_density %>% filter(year==2024)) + 
  geom_point(aes(x=hour, y=call_density), alpha=0.1) + 
  geom_smooth(aes(x=hour, y=call_density), span = 0.1) + 
  facet_wrap(~vegetation)

# Lesser Goldfinch
summary(lm(data = lego_density %>% 
             filter(year==2024), 
           call_density ~ sin((hour+3)*pi/24) * yday + vegetation))
ggplot(lego_density %>% filter(year==2024)) + 
  geom_point(aes(x=hour, y=call_density), alpha=0.1) + 
  geom_smooth(aes(x=hour, y=call_density), span = 0.1) + 
  facet_wrap(~vegetation)

# Spotted Towhee
summary(lm(data = spto_density %>% 
             filter(year==2024), 
           call_density ~ sin((hour+3)*pi/24) * yday + vegetation))
ggplot(spto_density %>% filter(year==2024)) + 
  geom_point(aes(x=hour, y=call_density), alpha=0.1) + 
  geom_smooth(aes(x=hour, y=call_density), span = 0.1) + 
  facet_wrap(~vegetation)

# Hooded Oriole
summary(lm(data = hoor_density %>% 
             filter(year==2024), 
           call_density ~ sin((hour+3)*pi/24) * yday + vegetation))
ggplot(hoor_density %>% filter(year==2024)) + 
  geom_point(aes(x=hour, y=call_density), alpha=0.1) + 
  geom_smooth(aes(x=hour, y=call_density), span = 0.1) + 
  facet_wrap(~vegetation)

# Anna's Hummingbird
summary(lm(data = anhu_density %>% 
             filter(year==2024), 
           call_density ~ sin((hour+3)*pi/24) * yday + vegetation))
ggplot(anhu_density %>% filter(year==2024)) + 
  geom_point(aes(x=hour, y=call_density), alpha=0.1) + 
  geom_smooth(aes(x=hour, y=call_density), span = 0.1) + 
  facet_wrap(~vegetation)

# Allen's Hummingbird
summary(lm(data = alhu_density %>% 
             filter(year==2024), 
           call_density ~ sin((hour+3)*pi/24) * yday + vegetation))
ggplot(alhu_density %>% filter(year==2024)) + 
  geom_point(aes(x=hour, y=call_density), alpha=0.1) + 
  geom_smooth(aes(x=hour, y=call_density), span = 0.1) + 
  facet_wrap(~vegetation)

# Build density data for each species in the guild list
species_density <- lapply(species_guilds$common_name,
                        getHourlyDensity)
species_density <- merge(bind_rows(species_density),
                       species_guilds,
                       by=("common_name"))
# Get vocal density by species guild over time / location 
guild_density <- species_density %>%
  group_by(site, year, yday, hour) %>%
  summarize(aerial_density = sum(call_density*aerial_insectivore),
            gleaning_density = sum(call_density*gleaning_insectivore),
            granivore_density = sum(call_density*granivore),
            carnivore_density = sum(call_density*carnivore)) %>%
  mutate(total_density = aerial_density + gleaning_density + granivore_density + carnivore_density)


# Build model comparing gleaning density by creek
gleaning_model <- lm(data = guild_density %>% 
                        filter(year==2024, 
                               substr(site,1,2) %in% c("AC","EW")), 
                      gleaning_density ~ sin((hour+3)*pi/24) * yday + (substr(site,1,2)=="AC"))
summary(gleaning_model)
predict(gleaning_model, 
        newdata = data.frame(yday = lubridate::yday(as.Date("2024-05-25")),
                             hour = 9, 
                             site = "AC_05"))
predict(gleaning_model, 
        newdata = data.frame(yday = lubridate::yday(as.Date("2024-05-25")),
                             hour = 9, 
                             site = "EW_01"))
# Visualize differences by time of day, season, and creek
ggplot(guild_density %>% filter(substr(site,1,2) %in% c("AC", "EW"), year==2024)) + 
  geom_point(aes(x=hour, y=gleaning_density), alpha=0.1) + 
  geom_smooth(aes(x=hour, y=gleaning_density), span = 0.1) + 
  facet_wrap(~substr(site,1,2))

# Build model comparing granivore density by creek
granivore_model <- lm(data = guild_density %>% 
                        filter(year==2024, 
                               substr(site,1,2) %in% c("AC","EW")), 
                      granivore_density ~ sin((hour+3)*pi/24) * yday + (substr(site,1,2)=="AC"))
summary(granivore_model)
predict(granivore_model, 
        newdata = data.frame(yday = lubridate::yday(as.Date("2024-05-25")),
                             hour = 9, 
                             site = "AC_05"))
predict(granivore_model, 
        newdata = data.frame(yday = lubridate::yday(as.Date("2024-05-25")),
                             hour = 9, 
                             site = "EW_01"))
# Visualize differences by time of day, season, and creek
ggplot(guild_density %>% filter(substr(site,1,2) %in% c("AC", "EW"), year==2024)) + 
  geom_point(aes(x=hour, y=granivore_density), alpha=0.1) + 
  geom_smooth(aes(x=hour, y=granivore_density), span = 0.1) + 
  facet_wrap(~substr(site,1,2))


