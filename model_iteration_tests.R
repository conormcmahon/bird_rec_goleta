
# Build several occupancy models each for multiple species
# Iterate parameters used between models
# Compare results and choose best fit models on basis of AICc
# Generate summary figures
# For now, using single-season occupancy models, building separate models for 2022 and 2023

# Suppress some warnings from dplyr's summarize() function
options(dplyr.summarise.inform = FALSE)

source(here::here("simple_presence_assessment.R"))

generateModels <- function(target_species, years, hour_offset, quiet=FALSE)
{
  
  # ****************************************************
  # ******************* Build Models ******************* 
  
  # ***************** 5 State Parameters ***************** 
  # Model with all considered variables
  first_fit <- modelSpecies(target_species,
                           ~doy + hour + surface_flow
                           ~wet + decid + height_pct_80 + low_veg + elevation,
                           rep(-1, 10), 
                           hour_offset = hour_offset,
                           use_point_counts = TRUE,
                           years_included = years,
                           quiet=TRUE)
  full_mod <- first_fit[[1]]
  known_sites = first_fit[[7]]
  covariate_data <- first_fit[[2]]
  
  # ***************** 4 State Parameters ***************** 
  wdhl_mod <- occu(data = covariate_data,
                   formula = ~doy + hour + surface_flow
                             ~wet + decid + height_pct_80 + low_veg,
                   starts = rep(-1, 9),
                   knownOcc = known_sites)
  
  wdhe_mod <- occu(data = covariate_data,
                   formula = ~doy + hour + surface_flow
                   ~wet + decid + height_pct_80 + elevation,
                   starts = rep(-1, 9),
                   knownOcc = known_sites)
  
  wdle_mod <- occu(data = covariate_data,
                   formula = ~doy + hour + surface_flow
                   ~wet + decid + low_veg + elevation,
                   starts = rep(-1, 9),
                   knownOcc = known_sites)
  
  wleh_mod <- occu(data = covariate_data,
                   formula = ~doy + hour + surface_flow
                   ~wet + low_veg + height_pct_80 + elevation,
                   starts = rep(-1, 9),
                   knownOcc = known_sites)
  
  dlhe_mod <- occu(data = covariate_data,
                   formula = ~doy + hour + surface_flow
                   ~ decid + low_veg + height_pct_80 + elevation,
                   starts = rep(-1, 9),
                   knownOcc = known_sites)

  # ***************** 3 State Parameters *****************
  
  wdl_mod <- occu(data = covariate_data,
                  formula = ~doy + hour + surface_flow
                  ~ wet + decid + low_veg,
                  starts = rep(-1, 8),
                  knownOcc = known_sites)

  wdh_mod <- occu(data = covariate_data,
                  formula = ~doy + hour + surface_flow
                  ~ wet + decid + height_pct_80,
                  starts = rep(-1, 8),
                  knownOcc = known_sites)

  wde_mod <- occu(data = covariate_data,
                  formula = ~doy + hour + surface_flow
                  ~ wet + decid + elevation,
                  starts = rep(-1, 8),
                  knownOcc = known_sites)

  wlh_mod <- occu(data = covariate_data,
                 formula = ~doy + hour + surface_flow
                 ~ wet + low_veg + height_pct_80,
                 starts = rep(-1, 8),
                 knownOcc = known_sites)

  wle_mod <- occu(data = covariate_data,
                  formula = ~doy + hour + surface_flow
                  ~ wet + low_veg + elevation,
                  starts = rep(-1, 8),
                  knownOcc = known_sites)

  whe_mod <- occu(data = covariate_data,
                  formula = ~doy + hour + surface_flow
                  ~ wet + height_pct_80 + elevation,
                  starts = rep(-1, 8),
                  knownOcc = known_sites)

  dlh_mod <- occu(data = covariate_data,
                  formula = ~doy + hour + surface_flow
                  ~ decid + low_veg + height_pct_80,
                  starts = rep(-1, 8),
                  knownOcc = known_sites)

  dle_mod <- occu(data = covariate_data,
                  formula = ~doy + hour + surface_flow
                  ~ decid + low_veg + elevation,
                  starts = rep(-1, 8),
                  knownOcc = known_sites)
  
  dhe_mod <- occu(data = covariate_data,
                  formula = ~doy + hour + surface_flow
                  ~ decid + height_pct_80 + elevation,
                  starts = rep(-1, 8),
                  knownOcc = known_sites)
  
  lhe_mod <- occu(data = covariate_data,
                  formula = ~doy + hour + surface_flow
                  ~ low_veg + height_pct_80 + elevation,
                  starts = rep(-1, 8),
                  knownOcc = known_sites)

  # ***************** 2 State Parameters *****************
  wd_mod <- occu(data = covariate_data,
                  formula = ~doy + hour + surface_flow
                  ~ wet + decid,
                  starts = rep(-1, 7),
                  knownOcc = known_sites)
  
  wl_mod <- occu(data = covariate_data,
                 formula = ~doy + hour + surface_flow
                 ~ wet + low_veg,
                 starts = rep(-1, 7),
                 knownOcc = known_sites)
  
  wh_mod <- occu(data = covariate_data,
                 formula = ~doy + hour + surface_flow
                 ~ wet + height_pct_80,
                 starts = rep(-1, 7),
                 knownOcc = known_sites)
  
  we_mod <- occu(data = covariate_data,
                 formula = ~doy + hour + surface_flow
                 ~ wet + elevation,
                 starts = rep(-1, 7),
                 knownOcc = known_sites)
  
  dl_mod <- occu(data = covariate_data,
                 formula = ~doy + hour + surface_flow
                 ~ decid + low_veg,
                 starts = rep(-1, 7),
                 knownOcc = known_sites)
  
  dh_mod <- occu(data = covariate_data,
                 formula = ~doy + hour + surface_flow
                 ~ decid + height_pct_80,
                 starts = rep(-1, 7),
                 knownOcc = known_sites)
  
  de_mod <- occu(data = covariate_data,
                 formula = ~doy + hour + surface_flow
                 ~ decid + elevation,
                 starts = rep(-1, 7),
                 knownOcc = known_sites)
  
  lh_mod <- occu(data = covariate_data,
                 formula = ~doy + hour + surface_flow
                 ~ low_veg + height_pct_80,
                 starts = rep(-1, 7),
                 knownOcc = known_sites)
  
  le_mod <- occu(data = covariate_data,
                 formula = ~doy + hour + surface_flow
                 ~ low_veg + elevation,
                 starts = rep(-1, 7),
                 knownOcc = known_sites)
  
  he_mod <- occu(data = covariate_data,
                 formula = ~doy + hour + surface_flow
                 ~ height_pct_80 + elevation,
                 starts = rep(-1, 7),
                 knownOcc = known_sites)
  
  # ***************** 2 State Parameters *****************
  
  w_mod <- occu(data = covariate_data,
                 formula = ~doy + hour + surface_flow
                 ~ wet,
                 starts = rep(-1, 6),
                 knownOcc = known_sites)
  
  d_mod <- occu(data = covariate_data,
                formula = ~doy + hour + surface_flow
                ~ decid,
                starts = rep(-1, 6),
                knownOcc = known_sites)
  
  l_mod <- occu(data = covariate_data,
                formula = ~doy + hour + surface_flow
                ~ low_veg,
                starts = rep(-1, 6),
                knownOcc = known_sites)
  
  h_mod <- occu(data = covariate_data,
                formula = ~doy + hour + surface_flow
                ~ height_pct_80,
                starts = rep(-1, 6),
                knownOcc = known_sites)
  
  e_mod <- occu(data = covariate_data,
                formula = ~doy + hour + surface_flow
                ~ elevation,
                starts = rep(-1, 6),
                knownOcc = known_sites)

  # ***************** 0 State Parameters *****************
  # Model with only intercepts for state formula, but still includes detection parameters
  det_mod <- occu(data = covariate_data,
                  formula = ~doy + hour + surface_flow
                  ~1,
                  starts = rep(-1, 5),
                  knownOcc = known_sites)

  # ***************** Full Null Model Parameters *****************
  # Model with only intercepts, even for detection
  null_mod <- occu(data = covariate_data,
                   formula = ~1
                   ~1,
                   starts = rep(-1, 2),
                   knownOcc = known_sites)

  # ****************************************************
  # ***************** Assemble Outputs *****************
  
  # Construct list with all models contained in it
  model_list <- list(full_mod = full_mod,
                     wdhl_mod = wdhl_mod,
                     wdhe_mod = wdhe_mod,
                     wdle_mod = wdle_mod,
                     wleh_mod = wleh_mod,
                     dlhe_mod = dlhe_mod,
                     wdl_mod = wdl_mod,
                     wdh_mod = wdh_mod,
                     wde_mod = wde_mod,
                     wlh_mod = wlh_mod,
                     wle_mod = wle_mod,
                     whe_mod = whe_mod,
                     dlh_mod = dlh_mod,
                     dle_mod = dle_mod,
                     dhe_mod = dhe_mod,
                     lhe_mod = lhe_mod,
                     wd_mod = wd_mod,
                     wl_mod = wl_mod,
                     wh_mod = wh_mod,
                     we_mod = we_mod,
                     dl_mod = dl_mod,
                     dh_mod = dh_mod,
                     de_mod = de_mod,
                     lh_mod = lh_mod,
                     le_mod = le_mod,
                     he_mod = he_mod,
                     w_mod = w_mod,
                     d_mod = d_mod,
                     l_mod = l_mod,
                     h_mod = h_mod,
                     e_mod = e_mod,
                     det_mod = det_mod,
                     null_mod = null_mod)
  
  model_comparison <- aictab(model_list) # remove extra data, retain only model
  
  # State coefficients for best model
  #   Coefficient estimates
  #   95% confidence intervals
  #   p-values and statistics
  best_state_coefs <- summary(model_list[[model_comparison$Modnames[1]]])$state
  conf_state <- confint(model_list[[model_comparison$Modnames[1]]], type="state")
  best_state_coefs <- cbind(best_state_coefs, conf_state)
  # Detection coefficients for best model
  #   Coefficient estimates
  #   95% confidence intervals
  #   p-values and statistics
  best_detection_coefs <- summary(model_list[[model_comparison$Modnames[1]]])$det
  conf_detection <- confint(model_list[[model_comparison$Modnames[1]]], type="det")
  best_detection_coefs <- cbind(best_detection_coefs, conf_detection)
  
  if(!quiet)
  {
    print(head(model_comparison))
    print(model_list[[model_comparison$Modnames[1]]])
  }
  
  return(list(model_comparison = model_comparison,
              all_models = model_list,
              state_coefs_best = best_state_coefs,
              detection_coefs_best = best_detection_coefs,
              data_input = covariate_data,
              species = target_species, 
              year = years))
}

# NOTE - did some manual iteration for each species to determine best hour offset...
# Iteratively ran models with only hour + doy for detection in 2022, slightly changing hour_offset
#   to select a value which minimized AIC
# In future, could come back to automate this, or even fit a better curve than sinusoid
getBestHour <- function(target_species, year)
{
  model_3 <- modelSpecies(target_species,
                          ~doy + hour
                          ~1,
                          rep(-1, 4), 
                          hour_offset = 3,
                          use_point_counts = TRUE,
                          years_included = year,
                          quiet=TRUE)
  model_2p5 <- modelSpecies(target_species,
                          ~doy + hour
                          ~1,
                          rep(-1, 4), 
                          hour_offset = 2.5,
                          use_point_counts = TRUE,
                          years_included = year,
                          quiet=TRUE)
  model_2 <- modelSpecies(target_species,
                          ~doy + hour
                          ~1,
                          rep(-1, 4), 
                          hour_offset = 2,
                          use_point_counts = TRUE,
                          years_included = year,
                          quiet=TRUE)
  model_1p5 <- modelSpecies(target_species,
                          ~doy + hour
                          ~1,
                          rep(-1, 4), 
                          hour_offset = 1.5,
                          use_point_counts = TRUE,
                          years_included = year,
                          quiet=TRUE)
  model_1 <- modelSpecies(target_species,
                          ~doy + hour
                          ~1,
                          rep(-1, 4), 
                          hour_offset = 1,
                          use_point_counts = TRUE,
                          years_included = year,
                          quiet=TRUE)
  model_p5 <- modelSpecies(target_species,
                          ~doy + hour
                          ~1,
                          rep(-1, 4), 
                          hour_offset = 0.5,
                          use_point_counts = TRUE,
                          years_included = year,
                          quiet=TRUE)
  model_0 <- modelSpecies(target_species,
                          ~doy + hour
                          ~1,
                          rep(-1, 4), 
                          hour_offset = 0,
                          use_point_counts = TRUE,
                          years_included = year,
                          quiet=TRUE)
  model_np5 <- modelSpecies(target_species,
                          ~doy + hour
                          ~1,
                          rep(-1, 4), 
                          hour_offset = -0.5,
                          use_point_counts = TRUE,
                          years_included = year,
                          quiet=TRUE)
  model_n1 <- modelSpecies(target_species,
                          ~doy + hour
                          ~1,
                          rep(-1, 4), 
                          hour_offset = -1,
                          use_point_counts = TRUE,
                          years_included = year,
                          quiet=TRUE)
  
  aic_results <- aictab(list(model_3=model_3[[1]],
                             model_2p5=model_2p5[[1]],
                             model_2=model_2[[1]],
                             model_1p5=model_1p5[[1]],
                             model_1=model_1[[1]],
                             model_p5=model_p5[[1]],
                             model_0=model_0[[1]],
                             model_np5=model_np5[[1]],
                             model_n1=model_n1[[1]]))
  hour_options <- (6:-2)/2
  names(hour_options) <- c("model_3",
                           "model_2p5",
                           "model_2",
                           "model_1p5",
                           "model_1",
                           "model_p5",
                           "model_0",
                           "model_np5",
                           "model_n1")
  return(list(best_hour_offset = as.numeric(hour_options[aic_results$Modnames[1]]),
              aic_results = aic_results,
              species = target_species))
}


# Black-headed Grosbeak
bhgr_hr_offset <- getBestHour("Black-headed Grosbeak", 2022)
bhgr_22 <- generateModels("Black-headed Grosbeak", 2022, bhgr_hr_offset[[1]])
bhgr_23 <- generateModels("Black-headed Grosbeak", 2023, bhgr_hr_offset[[1]])

# Song Sparrow
sosp_hr_offset <- getBestHour("Song Sparrow", 2022)
sosp_22 <- generateModels("Song Sparrow", 2022, sosp_hr_offset[[1]])
sosp_23 <- generateModels("Song Sparrow", 2023, sosp_hr_offset[[1]])

# Yellow Warbler
yewa_hr_offset <- getBestHour("Yellow Warbler", 2022)
yewa_22 <- generateModels("Yellow Warbler", 2022, yewa_hr_offset[[1]])
yewa_23 <- generateModels("Yellow Warbler", 2023, yewa_hr_offset[[1]])

# Wilson's Warbler
wiwa_hr_offset <- getBestHour("Wilson's Warbler", 2022)
wiwa_22 <- generateModels("Wilson's Warbler", 2022, wiwa_hr_offset[[1]])
wiwa_23 <- generateModels("Wilson's Warbler", 2023, wiwa_hr_offset[[1]])

# Warbling Vireo
wavi_hr_offset <- getBestHour("Warbling Vireo", 2022)
wavi_22 <- generateModels("Warbling Vireo", 2022, wavi_hr_offset[[1]])
wavi_23 <- generateModels("Warbling Vireo", 2023, wavi_hr_offset[[1]])

# Purple Finch
pufi_hr_offset <- getBestHour("Purple Finch", 2022)
pufi_22 <- generateModels("Purple Finch", 2022, pufi_hr_offset[[1]])
pufi_23 <- generateModels("Purple Finch", 2023, pufi_hr_offset[[1]])

# American Robin
amro_hr_offset <- getBestHour("American Robin", 2022)
amro_22 <- generateModels("American Robin", 2022, amro_hr_offset[[1]])
amro_23 <- generateModels("American Robin", 2023, amro_hr_offset[[1]])

# Chestnut-backed Chickadee
cbch_hr_offset <- getBestHour("Chestnut-backed Chickadee", 2022)
cbch_22 <- generateModels("Chestnut-backed Chickadee", 2022, cbch_hr_offset[[1]])
cbch_23 <- generateModels("Chestnut-backed Chickadee", 2023, cbch_hr_offset[[1]])

# Black Phoebe
blph_hr_offset <- getBestHour("Black Phoebe", 2022)
blph_22 <- generateModels("Black Phoebe", 2022, blph_hr_offset[[1]])
cbch_23 <- generateModels("Black Phoebe", 2023, blph_hr_offset[[1]])

# Common Yellowthroat
coye_hr_offset <- getBestHour("Common Yellowthroat", 2022)
coye_22 <- generateModels("Common Yellowthroat", 2022, coye_hr_offset[[1]])
coye_23 <- generateModels("Common Yellowthroat", 2023, coye_hr_offset[[1]])

# Canyon Wren
cawr_hr_offset <- getBestHour("Canyon Wren", 2022)
cawr_22 <- generateModels("Canyon Wren", 2022, cawr_hr_offset[[1]])
cawr_23 <- generateModels("Canyon Wren", 2023, cawr_hr_offset[[1]])

# Hooded Oriole
hoor_hr_offset <- getBestHour("Hooded Oriole", 2022)
hoor_22 <- generateModels("Hooded Oriole", 2022, hoor_hr_offset[[1]])
hoor_23 <- generateModels("Hooded Oriole", 2023, hoor_hr_offset[[1]])

# Blue Grosbeak
blgr_hr_offset <- getBestHour("Blue Grosbeak", 2022)
blgr_22 <- generateModels("Blue Grosbeak", 2022, blgr_hr_offset[[1]])
blgr_23 <- generateModels("Blue Grosbeak", 2023, blgr_hr_offset[[1]])

# Bullock's Oriole
buor_hr_offset <- getBestHour("Bullock's Oriole", 2022)
buor_22 <- generateModels("Bullock's Oriole", 2022, buor_hr_offset[[1]])
buor_23 <- generateModels("Bullock's Oriole", 2023, buor_hr_offset[[1]])

# Pacific-slope Flycatcher
psfl_hr_offset <- getBestHour("Pacific-slope Flycatcher", 2022)
psfl_22 <- generateModels("Pacific-slope Flycatcher", 2022, psfl_hr_offset[[1]])
psfl_23 <- generateModels("Pacific-slope Flycatcher", 2023, psfl_hr_offset[[1]])

# Spotted Towhee
spto_hr_offset <- getBestHour("Spotted Towhee", 2022)
spto_22 <- generateModels("Spotted Towhee", 2022, spto_hr_offset[[1]])
spto_23 <- generateModels("Spotted Towhee", 2023, spto_hr_offset[[1]])

# Oak Titmouse
oati_hr_offset <- getBestHour("Oak Titmouse", 2022)
oati_22 <- generateModels("Oak Titmouse", 2022, oati_hr_offset[[1]])
oati_23 <- generateModels("Oak Titmouse", 2023, oati_hr_offset[[1]])

# Acorn Woodpecker
acwo_hr_offset <- getBestHour("Acorn Woodpecker", 2022)
acwo_22 <- generateModels("Acorn Woodpecker", 2022, acwo_hr_offset[[1]])
acwo_23 <- generateModels("Acorn Woodpecker", 2023, acwo_hr_offset[[1]])

# Orange-crowned Warbler
ocwa_hr_offset <- getBestHour("Orange-crowned Warbler", 2022)
ocwa_22 <- generateModels("Orange-crowned Warbler", 2022, ocwa_hr_offset[[1]])
ocwa_23 <- generateModels("Orange-crowned Warbler", 2023, ocwa_hr_offset[[1]])

# American Crow
amcr_hr_offset <- getBestHour("American Crow", 2022)
amcr_22 <- generateModels("American Crow", 2022, amcr_hr_offset[[1]])
amcr_23 <- generateModels("American Crow", 2023, amcr_hr_offset[[1]])

# House Finch
hofi_hr_offset <- getBestHour("House Finch", 2022)
hofi_22 <- generateModels("House Finch", 2022, hofi_hr_offset[[1]])
hofi_23 <- generateModels("House Finch", 2023, hofi_hr_offset[[1]])

# Lesser Goldfinch
lego_hr_offset <- getBestHour("Lesser Goldfinch", 2022)
lego_22 <- generateModels("Lesser Goldfinch", 2022, lego_hr_offset[[1]])
lego_23 <- generateModels("Lesser Goldfinch", 2023, lego_hr_offset[[1]])

# California Scrub-Jay
casj_hr_offset <- getBestHour("California Scrub-Jay", 2022)
casj_22 <- generateModels("California Scrub-Jay", 2022, casj_hr_offset[[1]])
casj_23 <- generateModels("California Scrub-Jay", 2023, casj_hr_offset[[1]])

# Wrentit
wren_hr_offset <- getBestHour("Wrentit", 2022)
wren_22 <- generateModels("Wrentit", 2022, wren_hr_offset[[1]])
wren_23 <- generateModels("Wrentit", 2023, wren_hr_offset[[1]])

# Anna's Hummingbird
anhu_hr_offset <- getBestHour("Anna's Hummingbird", 2022)
anhu_22 <- generateModels("Anna's Hummingbird", 2022, anhu_hr_offset[[1]])
anhu_23 <- generateModels("Anna's Hummingbird", 2023, anhu_hr_offset[[1]])

# Bewick's Wren
bewr_hr_offset <- getBestHour("Bewick's Wren", 2022)
bewr_22 <- generateModels("Bewick's Wren", 2022, bewr_hr_offset[[1]])
bewr_23 <- generateModels("Bewick's Wren", 2023, bewr_hr_offset[[1]])

# Mourning Dove
modo_hr_offset <- getBestHour("Mourning Dove", 2022)
modo_22 <- generateModels("Mourning Dove", 2022, modo_hr_offset[[1]])
modo_23 <- generateModels("Mourning Dove", 2023, modo_hr_offset[[1]])

# Bushtit
bush_hr_offset <- getBestHour("Bushtit", 2022)
bush_22 <- generateModels("Bushtit", 2022, bush_hr_offset[[1]])
bush_23 <- generateModels("Bushtit", 2023, bush_hr_offset[[1]])

# Dark-eyed Junco
deju_hr_offset <- getBestHour("Dark-eyed Junco", 2022)
deju_22 <- generateModels("Dark-eyed Junco", 2022, deju_hr_offset[[1]])
deju_23 <- generateModels("Dark-eyed Junco", 2023, deju_hr_offset[[1]])

# White-breasted Nuthatch
wbnu_hr_offset <- getBestHour("White-breasted Nuthatch", 2022)
wbnu_22 <- generateModels("White-breasted Nuthatch", 2022, wbnu_hr_offset[[1]])
wbnu_23 <- generateModels("White-breasted Nuthatch", 2023, wbnu_hr_offset[[1]])

# Western Bluebird
webl_hr_offset <- getBestHour("Western Bluebird", 2022)
webl_22 <- generateModels("Western Bluebird", 2022, webl_hr_offset[[1]])
webl_23 <- generateModels("Western Bluebird", 2023, webl_hr_offset[[1]])



# Several meta-tests to try...
#  1) Linear regression, decid and wet coefficients from 2022 vs. 2023 by species 
#      does wet coefficient weaken when water is everywhere?
#      compare values from either dw_mod, or dw_mod, OR from best fit and set coefficient to 0 when unused
#  2) Color-scaled table of coefficient values for best-fit models, with 0s for unused variables
#  3) Scatter plot + SE Error Bars for several species in 2022 and 2023, with x=coef_wet, y=coef_decid 
# Detection vs. Year / DOY
# Map of predicted habitat for YEWA, WIWA, WAVI? 
# Predicted abundance by species for both years? 


# Compare water and deciduous coefficients by species for 2022 and 2023
getComparison <- function(target_data, species=target_data$species, year=target_data$year, model_choice="wd_mod")
{
  target_model <- target_data[[2]][[model_choice]]
  coef_summary <- summary(target_model)
  
  # State coefficients and confidence intervals for best model
  state_coefs <- summary(target_model)$state
  conf_state <- confint(target_model, type="state")
  state_coefs <- cbind(state_coefs, conf_state)
  # Detection coefficients and confidence intervals for best model
  detection_coefs <- summary(target_model)$det
  conf_detection <- confint(target_model, type="det")
  detection_coefs <- cbind(detection_coefs, conf_detection)
  
  output_df <- data.frame(rbind(state_coefs,
                                detection_coefs))
  output_df$species <- species
  output_df$year <- year
  
  output_df$variable <- rownames(output_df)
  rownames(output_df) <- NULL
  
  return(output_df)
}


createGroupedPlots <- function(model_name, variable_name, label_prefix, plot_limits=c(-3,6))
{
  # Create list of all species models
  model_list <- list(bhgr_22, bhgr_23,
                     sosp_22, sosp_23,
                     wiwa_22, wiwa_23,
                     yewa_22, yewa_23,
                     pufi_22, pufi_23,
                     cbch_22, cbch_23,
                     spto_22, spto_23,
                     psfl_22, psfl_23,
                     oati_22, oati_23,
                     acwo_22, acwo_23,
                     amcr_22, amcr_23,
                     hofi_22, hofi_23,
                     lego_22, lego_23,
                     casj_22, casj_23,
                     wren_22, wren_23,
                     anhu_22, anhu_23,
                     amro_22, amro_23,
                     blph_22, blph_23, 
                     coye_22, coye_23, 
                     bewr_22, bewr_23, 
                     modo_22, modo_23, 
                     bush_22, bush_23, 
                     deju_22, deju_23, 
                     wbnu_22, wbnu_23, 
                     cawr_22, cawr_23, 
                     hoor_22, hoor_23, 
                     blgr_22, blgr_23,
                     buor_22, buor_23,
                     webl_22, webl_23)
  # Extract coefficient data for each model
  coef_comparison <- bind_rows(lapply(model_list,
                                      getComparison,
                                      model_choice=model_name))
  # Add year to dataframe
  coef_comparison$year <- rep(c(rep(2022,nrow(coef_comparison)/length(model_list)), 
                                rep(2023,nrow(coef_comparison)/length(model_list))), 
                              length(model_list)/2)
  # Add species to dataframe
  coef_comparison$species <- c(rep("Black-headed Grosbeak", nrow(coef_comparison)/length(model_list)*2),
                               rep("Song Sparrow", nrow(coef_comparison)/length(model_list)*2),
                               rep("Wilson's Warbler", nrow(coef_comparison)/length(model_list)*2),
                               rep("Yellow Warbler", nrow(coef_comparison)/length(model_list)*2),
                               rep("Purple Finch", nrow(coef_comparison)/length(model_list)*2),
                               rep("Chestnut-backed Chickadee", nrow(coef_comparison)/length(model_list)*2),
                               rep("Spotted Towhee", nrow(coef_comparison)/length(model_list)*2),
                               rep("Pacific-slope Flycatcher", nrow(coef_comparison)/length(model_list)*2),
                               rep("Oak Titmouse", nrow(coef_comparison)/length(model_list)*2),
                               rep("Acorn Woodpecker", nrow(coef_comparison)/length(model_list)*2),
                               rep("American Crow", nrow(coef_comparison)/length(model_list)*2),
                               rep("House Finch", nrow(coef_comparison)/length(model_list)*2),
                               rep("Lesser Goldfinch", nrow(coef_comparison)/length(model_list)*2),
                               rep("California Scrub-Jay", nrow(coef_comparison)/length(model_list)*2),
                               rep("Wrentit", nrow(coef_comparison)/length(model_list)*2),
                               rep("Anna's Hummingbird", nrow(coef_comparison)/length(model_list)*2),
                               rep("American Robin", nrow(coef_comparison)/length(model_list)*2),
                               rep("Black Phoebe", nrow(coef_comparison)/length(model_list)*2),
                               rep("Common Yellowthroat", nrow(coef_comparison)/length(model_list)*2),
                               rep("Bewick's Wren", nrow(coef_comparison)/length(model_list)*2),
                               rep("Mourning Dove", nrow(coef_comparison)/length(model_list)*2),
                               rep("Bushtit", nrow(coef_comparison)/length(model_list)*2),
                               rep("Dark-eyed Junco", nrow(coef_comparison)/length(model_list)*2),
                               rep("White-breasted Nuthatch", nrow(coef_comparison)/length(model_list)*2),
                               rep("Canyon Wren", nrow(coef_comparison)/length(model_list)*2),
                               rep("Hooded Oriole", nrow(coef_comparison)/length(model_list)*2),
                               rep("Blue Grosbeak", nrow(coef_comparison)/length(model_list)*2),
                               rep("Bullock's Oriole", nrow(coef_comparison)/length(model_list)*2),
                               rep("Western Bluebird", nrow(coef_comparison)/length(model_list)*2))
  
  # Get wide dataframe where rows are species, columns have coefficient estimate in 2022 and 2023
  wide_22_23 <- coef_comparison %>% 
    filter(variable %in% c(variable_name)) %>% 
    dplyr::select(Estimate, variable, year, species) %>% 
    pivot_wider(names_from=year, values_from=Estimate)
  names(wide_22_23) <- c("variable", "species", "v_2022", "v_2023")
  print(summary(lm(data=wide_22_23, v_2023 ~ v_2022)))
  
  # Get wide dataframe where rows are species, columns have standard error of estimate in 2022 and 2023
  wide_se <- coef_comparison %>% filter(variable %in% c(variable_name)) %>% dplyr::select(SE, variable, year, species) %>% pivot_wider(names_from=year, values_from=SE)
  names(wide_se) <- c("variable", "species", "v_2022", "v_2023")
  summary(lm(data=wide_se, v_2023 ~ v_2022))
  wide_se
  
  # Combine wide dataframes with coefficients and SE
  wide_22_23$se_2022 <- wide_se$v_2022
  wide_22_23$se_2023 <- wide_se$v_2023
  wide_22_23
  
  # List of obligately-riparian birds 
  riparian_birds <-  c("Black-headed Grosbeak", "Wilson's Warbler", "Warbling Vireo", "Yellow Warbler", "Song Sparrow", "Purple Finch",
                       "American Robin", "Common Yellowthroat", "Blue Grosbeak", "Chestnut-backed Chickadee")

  # Add boolean flag to data marking riparian species   
  wide_22_23$riparian <- wide_22_23$species %in% riparian_birds
  
  # Remove models which did not fit well 
  #   These might have very high standard error, or 
  #   NA values for SE or coefficients
  wide_data_filtered <- wide_22_23 %>% 
    drop_na(v_2022, v_2023, se_2022, se_2023) %>% 
    filter(se_2022 < 10, se_2023 < 10)
  
  # Get mean and SE of mean across filtered models 
  wide_data_stats <- wide_data_filtered %>% 
    group_by(riparian) %>% 
    summarize(mean_22 = mean(v_2022), 
              mean_23 = mean(v_2023),
              SE_22 = sqrt(sum(se_2022^2/n()^2)),
              SE_23 = sqrt(sum(se_2023^2/n()^2)))
  
  multiyear_plot <- ggplot() + 
    # All points
    geom_point(data=wide_data_filtered,
               aes(x=v_2022,
                   y=v_2023, 
                   col=(species %in% riparian_birds)), alpha=1) + 
    # Horizontal error bars for mean-of-models
    geom_pointrange(data=wide_data_stats,                       
                    aes(x=mean_22, y=mean_23, 
                        xmin=mean_22-SE_22, xmax=mean_22+SE_22, 
                        col=(riparian)),
                    size=1, linewidth=1) + 
    # Vertical error bars for mean-of-models
    geom_pointrange(data=wide_data_stats,
                    aes(x=mean_22, y=mean_23, 
                        ymin=mean_23-SE_23, ymax=mean_23+SE_23, 
                        col=(riparian)),
                    size=1, linewidth=1) + 
    # Horizontal and Vertical intercept lines
    geom_hline(yintercept = 0, col="black", linetype=1, alpha=0.2) + 
    geom_vline(xintercept = 0, col="black", linetype=1, alpha=0.2) + 
    # Theme and aesthetics
    theme_bw() + 
    scale_x_continuous(limits=plot_limits) + 
    scale_y_continuous(limits=plot_limits) + 
    theme(legend.position = "None") + 
    xlab(paste(label_prefix, " (2022)",sep="")) + 
    ylab(paste(label_prefix, " (2023)",sep=""))
  
  metamodel_comparison <- coef_comparison %>% 
    mutate(riparian = species %in% riparian_birds) %>%
    filter(variable == variable_name,
           #P...z.. < 0.1
           SE < abs(Estimate)) %>%
    group_by(riparian) %>% 
    summarize(mean_sens = mean(mean(Estimate)),
              std_err = sqrt(sum(SE^2/n()^2)),
              count = n())
  print(metamodel_comparison)
  
  return(list(plot = multiyear_plot,
              model_comparison = metamodel_comparison,
              all_data = coef_comparison))
}

# First, create plots relating bird presence to deciduous tree cover
grouped_decid_plot <- createGroupedPlots("full_mod", "decid", "Deciduous Trees", c(-5,10))
ggsave("D:/birdrec/reports/grouped_deciduous_response.pdf",
       grouped_decid_plot[[1]],
       width = 6.25, height = 6.25)
ggsave("D:/birdrec/reports/grouped_deciduous_response.png",
       grouped_decid_plot[[1]],
       width = 6.25, height = 6.25)

# Next, the same for surface water 
grouped_water_plot <- createGroupedPlots("full_mod", "wetTRUE", "Surface Water", c(-2.5,5))
ggsave("D:/birdrec/reports/grouped_water_response.pdf",
       grouped_water_plot[[1]],
       width = 6.25, height = 6.25)
ggsave("D:/birdrec/reports/grouped_water_response.png",
       grouped_water_plot[[1]],
       width = 6.25, height = 6.25)


# Get model summaries for all species
getCoefSummaries <- function(target_model_list, year, species)
{
  best_model_name <- target_model_list[[1]]$Modnames[1]
  best_model <- target_model_list[[2]][[best_model_name]]
  model_summary <- summary(best_model)
  summary_df <- as.data.frame(model_summary$state)
  names(summary_df) <- c("estimate", "se", "z", "p")
  summary_df <- summary_df %>% 
    mutate(coef_sig = estimate * (p < 0.1)) %>%
    dplyr::select(coef_sig) 
  summary_df <- data.frame(t(summary_df)) %>%
    mutate(year = year,
           species = species)
  return(summary_df)
}

# Generate best-fit model summaries for a bunch of species of interest
yewa_22_coefs <- getCoefSummaries(yewa_22, 2022, "Yellow Warbler")
wiwa_22_coefs <- getCoefSummaries(wiwa_22, 2022, "Wilson's Warbler")
bhgr_22_coefs <- getCoefSummaries(bhgr_22, 2022, "Black-headed Grosbeak")
sosp_22_coefs <- getCoefSummaries(sosp_22, 2022, "Song Sparrow")
cbch_22_coefs <- getCoefSummaries(cbch_22, 2022, "Chestnut-backed Chickadee")
spto_22_coefs <- getCoefSummaries(spto_22, 2022, "Spotted Towhee")
acwo_22_coefs <- getCoefSummaries(acwo_22, 2022, "Acorn Woodpecker")
amcr_22_coefs <- getCoefSummaries(amcr_22, 2022, "American Crow")
ocwa_22_coefs <- getCoefSummaries(ocwa_22, 2022, "Orange-crowned Warbler")
psfl_22_coefs <- getCoefSummaries(psfl_22, 2022, "Western Flycatcher")
hofi_22_coefs <- getCoefSummaries(hofi_22, 2022, "House Finch")
modo_22_coefs <- getCoefSummaries(modo_22, 2022, "Mourning Dove")
wren_22_coefs <- getCoefSummaries(wren_22, 2022, "Wrentit")
blph_22_coefs <- getCoefSummaries(blph_22, 2022, "Black Phoebe")
deju_22_coefs <- getCoefSummaries(deju_22, 2022, "Dark-eyed Junco")
wbnu_22_coefs <- getCoefSummaries(wbnu_22, 2022, "White-breasted Nuthatch")
hoor_22_coefs <- getCoefSummaries(hoor_22, 2022, "Hooded Oriole")

# Mash all those dataframes into one big model summary table
all_coefs <- merge(yewa_22_coefs,
                   merge(wiwa_22_coefs,
                         merge(bhgr_22_coefs,
                               merge(sosp_22_coefs,
                                     merge(cbch_22_coefs,
                                           merge(spto_22_coefs,
                                                 merge(acwo_22_coefs,
                                                       merge(ocwa_22_coefs,
                                                             merge(psfl_22_coefs,
                                                                   merge(hofi_22_coefs,
                                                                         merge(modo_22_coefs,
                                                                               merge(wren_22_coefs,
                                                                                     merge(blph_22_coefs,
                                                                                           merge(deju_22_coefs,
                                                                                                 merge(wbnu_22_coefs,
                                                                                                       hoor_22_coefs, all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE),
                                                                   all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE), all=TRUE) %>%
  dplyr::select(species, year, X.Intercept., wetTRUE, decid, low_veg, height_pct_80, elevation)
all_coefs

# Show change in predicted occupancy between years
#   For most species, this is a small change - prob don't include in paper / poster
wiwa_22_pred <- predict(wiwa_22[[2]][[wiwa_22[[1]]$Modnames[1]]], newdata=wiwa_22$data_input, type='state')
wiwa_23_pred <- predict(wiwa_23[[2]][[wiwa_23[[1]]$Modnames[1]]], newdata=wiwa_23$data_input, type='state')
ggplot(data.frame(predicted = c(wiwa_22_pred$Predicted, wiwa_23_pred$Predicted),
                  year = c(rep(2022,112),rep(2023,103)))) + 
  geom_freqpoly(aes(x=predicted, group=year, col=year))

# Examples - predict presence and detection by vegetation, etc
#   Build simulation data for state variable
state_df <- data.frame(wet = c(rep(0,50),rep(1,50)),
                       decid = rep((0:49)/49,2),
                       low_veg = rep(0.6,100),
                       height_pct_80 = rep(0.7,100), 
                       elevation = rep(0.0,100))
#   Predict occupancy for several species
occupancy_prediction <- predict(bhgr_22[[2]][["full_mod"]], # model
                                newdata = state_df,
                                type="state")
state_df <- cbind(state_df, occupancy_prediction)
# Plot Occupancy
occupancy_plot <- ggplot(state_df) + 
  geom_line(aes(x=decid, y=Predicted, group=wet, col=factor(wet)), size=1) + 
  geom_line(aes(x=decid, y=Predicted-SE, group=wet, col=factor(wet)), linetype=2, size=0.8) + 
  geom_line(aes(x=decid, y=Predicted+SE, group=wet, col=factor(wet)), linetype=2, size=0.8) + 
  theme_bw() + 
  scale_x_continuous(expand=c(0,0), limits=c(0,1)) + 
  scale_y_continuous(expand=c(0,0), limits=c(0,1))
occupancy_plot
ggsave("D:/birdrec/reports/bhgr_occupancy_plot.png",
       occupancy_plot,
       width=8,height=6)

#   Now, the same for detection module. Start with data frame:
detection_df <- data.frame(doy = c(rep(0,50),rep(1,50)),
                       hour_real = rep((0:49)/49*24,2),
                       surface_flow = 0) %>%
  mutate(hour = sin(((hour_real+1)%%24)/24*pi),
         date = as.Date("2022-01-01") + doy*(223-164)+165)

#   Now, the same for detection module. Start with data frame:
detection_df <- data.frame(doy = c(rep(0,50),rep(1,50)),
                           hour_real = rep((0:49)/49*24,2),
                           surface_flow = 0) %>%
  mutate(hour = sin(((hour_real+1)%%24)/24*pi),
         date = as.Date("2022-01-01") + doy*(223-164)+165)

#   Predict occupancy for several species
detection_prediction <- predict(bhgr_22[[2]][["full_mod"]], # model
                                newdata = detection_df,
                                type="det")
detection_df <- cbind(detection_df, detection_prediction)
# Plot Occupancy
detection_plot <- ggplot(detection_df) + 
  geom_line(aes(x=hour_real, y=Predicted, group=doy, col=factor(doy)), size=1) + 
  geom_line(aes(x=hour_real, y=Predicted-SE, group=doy, col=factor(doy)), linetype=2, size=0.8) + 
  geom_line(aes(x=hour_real, y=Predicted+SE, group=doy, col=factor(doy)), linetype=2, size=0.8) +
  theme_bw() + 
  scale_x_continuous(expand=c(0,0), limits=c(0,24)) + 
  scale_y_continuous(expand=c(0,0), limits=c(0,1)) + 
  scale_color_manual(values = c("skyblue3", "gray50"),
                     labels = c(0, 1))
detection_plot
ggsave("D:/birdrec/reports/bhgr_detection_plot.png",
       detection_plot,
       width=8,height=6)

