
library(tidyverse)



# Surface water!
field_data <- st_read("D:/birdrec/field_data_goleta_2022/goleta_qfield_final_2023/birdrec_site_goleta_2022.gpkg") %>%
  drop_na(timestamp_collected) %>%                             # remove stolen or missing data 
  mutate(year = as.numeric(substr(timestamp_collected, 1, 4)), # add year
         site = str_replace(site_name, "-", "_")) %>%          # clean up formatting
  dplyr::select(-site_name)                                    # clean up formatting 
# Replace NA with FALSE
field_data$surface_flow <- replace_na(field_data$surface_flow, FALSE)
field_data$pooling_water <- replace_na(field_data$pooling_water, FALSE)
field_data$seeps_or_springs <- replace_na(field_data$seeps_or_springs, FALSE)
# Combine types of water into one surface water class
field_data$wet <- (field_data$surface_flow + field_data$pooling_water + field_data$seeps_or_springs) > 0
# Get variable to show whether water is ALWAYS present
field_data <- merge(as.data.frame(field_data), 
                    as.data.frame(field_data) %>% 
                      dplyr::select(-geometry) %>% 
                      group_by(site) %>% 
                      summarize(perennially_wet = sum(wet)==n()))
surface_water <- as.data.frame(field_data) %>% 
  dplyr::select(site, year, wet, perennially_wet) %>% 
  group_by(site, year) %>% 
  summarize(wet = sum(wet)>0,
            perennially_wet = sum(perennially_wet)==n())



# Remote Sensing
#   Loading this data from where it was constructed in ./remote_sensing_sampler.R
#   Come back and clean this up a lot later!
rs_data <- read_csv("D:/birdrec/remote_sensing_data/site_rs_data.csv") %>%
  mutate(site = site_name) %>%
  dplyr::select(-site_name)


# Field Vegetation Surveys
field_canopy_cover <- read_csv("D:/birdrec/vegetation/canopy_cover.csv") %>% 
  drop_na(Location, Site) %>%
  janitor::clean_names()
creeknames <- c("CSC",
                "SR",
                "NCOS",
                "KR",
                "MC",
                "AC",
                "SAC",
                "SJC",
                "AC",
                "EW",
                "LLC",
                "MC",
                "MC",
                "MC",
                "RC",
                "BR")
names(creeknames) <- c("Cold Spring Creek",
                       "Sedgwick Reserve",
                       "NCOS", 
                       "Kinevan Road", 
                       "Mission Creek (Botanic Garden)",
                       "Atascadero Creek",
                       "San Antonio Creek",
                       "San Jose Creek",
                       "Maria Ygnacio Creek (labeled AC)", 
                       "Ellwood",
                       "Lake Los Carneros",
                       "Mission Creek",
                       "Mission Creek (Rocky Nook)",
                       "Mission/Rattlesnake Confluence",
                       "Rattlesnake Creek",
                       "Baron Ranch")
field_canopy_cover$site <- paste(creeknames[field_canopy_cover$location], 
                                 sprintf("%02d", field_canopy_cover$site),
                                 sep="_")
field_canopy_summary <- field_canopy_cover %>%
  group_by(site, date) %>%
  summarize(willows = sum(total_cover_percent*(genus=="Salix")),
            deciduous_tree = (sum(plane_cover_percent*(genus %in% c("Platanus","Salix","Populus","Fraxinus","Acer","Alnus","Juglans","PLatanus"))) + 
              sum(plane_cover_percent*(species == "salicifolia")))/100,
            native_evergreen_tree = (sum(total_cover_percent*(genus %in% c("Umbellularia", "Pinus","Sequoia","Calocedrus","Notholithocarpus"))) + 
              sum(total_cover_percent*(species == "agrifolia")))/100,
            nonnative_evergreen_tree = sum(total_cover_percent*(genus %in% c("Acacia","Bauhinia","Eucalyptus","Evergreen broadleaf","Ficus","Jacaranda","Ligustrum","Nicotiana","Palm","Persea","Phoenix","Pittosporum","Prunus","Schinus")))/100,
            fruit_tree = sum(total_cover_percent*(genus %in% c("Ficus","Frangula","Heteromeles","Ligustrum","Manzanita","Melosma","Palm","Persea","Phoenoix","Pittosporum","Prunus","Rubus","Sambucus","Schinus","Toxicodendron")))/100,
            tall_tree_cover = sum(plane_cover_percent*(genus %in% c("Platanus", "Quercus","Pinus","Populus","Eucalyptus")))/100) %>%
  group_by(site) %>% 
  summarize(willows = mean(willows, na.rm=TRUE),
            deciduous_tree = mean(deciduous_tree, na.rm=TRUE),
            native_evergreen_tree = mean(native_evergreen_tree, na.rm=TRUE),
            nonnative_evergreen_tree = mean(nonnative_evergreen_tree, na.rm=TRUE),
            fruit_tree = mean(fruit_tree, na.rm=TRUE),
            tall_tree_cover = mean(tall_tree_cover, na.rm=TRUE))
# NOTE - we actually only have veg surveys at ~57 / 110 sites, so we won't use them in the main analysis, but we can use them to validate the remote sensing
rs_field_comparison <- merge(rs_data, field_canopy_summary)

# Compare field deciduous metric to remote sensing
decid_rs_model <- (lm(rs_field_comparison$deciduous_tree ~ rs_field_comparison$decid))
print(summary(decid_rs_field_model))
deciduous_plot <- ggplot() + 
  geom_point(data=rs_field_comparison,
             aes(x=decid, y=deciduous_tree)) + 
  geom_abline(slope = summary(decid_rs_model)$coef[2,1],
              intercept = summary(decid_rs_model)$coef[1,1]) + 
  geom_abline(slope = 1,
              intercept = 0, 
              col="black", linetype=2) + 
  scale_x_continuous(expand=c(0,0), limits=c(0,1)) + 
  scale_y_continuous(expand=c(0,0), limits=c(0,1)) + 
  theme_bw() +
  xlab("Deciduous Cover - Satellite Estimate") + 
  ylab("Deciduous Cover - Field Estimate")
deciduous_plot
ggsave("D:/birdrec/reports/rs_deciduousness_plot.png",
       deciduous_plot,
       width=6, height=6)

# Compare field estimates of tall trees to remote sensing estimate
tall_trees_rs_model <- (lm(rs_field_comparison$tall_tree_cover ~ rs_field_comparison$height_pct_80))
print(summary(tall_trees_rs_model))
tall_plot <- ggplot() + 
  geom_point(data=rs_field_comparison,
             aes(x=height_pct_80/max(height_pct_80), y=tall_tree_cover)) + 
  geom_abline(slope = summary(tall_trees_rs_model)$coef[2,1]*max(rs_field_comparison$height_pct_80),
              intercept = summary(tall_trees_rs_model)$coef[1,1]) +
  geom_abline(slope = 1,
              intercept = 0, 
              col="black", linetype=2) +
  theme_bw() + 
  xlab("80th Percentile of Tree Height") + 
  scale_x_continuous(expand=c(0,0), limits=c(0,1)) + 
  scale_y_continuous(expand=c(0,0), limits=c(0,1)) + 
  ylab("Fractional Cover by Tall Trees")
tall_plot
ggsave("D:/birdrec/reports/rs_tree_cover_plot.png",
       tall_plot,
       width=6, height=6) 





rs_water_comparison <- merge(rs_data, field_data)

# Deciduous - not predicted well by elevation or water 
ggplot(rs_water_comparison) + 
  geom_point(aes(x=wet, y=decid))
ggplot(rs_water_comparison) + 
  geom_point(aes(x=elevation, y=decid))
summary(lm(decid ~ wet+elevation, data=rs_water_comparison))
# LiDAR metrics
# Max height is positively correlated with std_max_height (R^2=0.41, coef = 0.179)
ggplot(rs_water_comparison) + 
  geom_point(aes(x=max_height, y=std_max_height))
summary(lm(std_max_height ~ max_height, data=rs_water_comparison))
# Max height is also, unsurprisingly, correlated with height percentile values
ggplot(rs_water_comparison) + 
  geom_point(aes(x=max_height, y=std_max_height))
summary(lm(std_max_height ~ max_height, data=rs_water_comparison))





