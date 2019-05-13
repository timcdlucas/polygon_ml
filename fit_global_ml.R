########
# master script for point Vs polygon Vs joint analysis
# Tim Lucas
# 2018-05-30
###########

if(Sys.info()["user"] != 'anita'){
  setwd('~/timz/timothy/polygon_ml')
} else {
  #setwd('~/Z/users/anita/point_polygon_join_comparison_analysis')
}

source("setUserInfo.R")

# define paths

PR_path <- Z('GBD2017/Processing/Stages/04b_PR_DB_Import_Export/Verified_Outputs/2018_02_15/pfpr_dump.csv')
API_path <- Z('GBD2017/Processing/Stages/04c_API_Data_Export/Checkpoint_Outputs/subnational.csv')
pop_path <- Z('GBD2017/Processing/Stages/03_Muster_Population_Figures/Verified_Outputs/Output_Pop_At_Risk_Pf_5K/ihme_corrected_frankenpop_All_Ages_3_2015_at_risk_pf.tif')
shapefile_path <- Z('master_geometries/Admin_Units/Global/GBD/GBD2017_MAP/GBD2017_MAP_MG_5K/')


cov_raster_paths <- c(
  Z('mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/5km/Synoptic/LST_Day_v6.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/EVI_v6/5km/Synoptic/EVI_v6.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Synoptic/TSI-Martens2-Pf.Synoptic.Overall.Mean.5km.Data.tif'),
  Z('GBD2017/Processing/Static_Covariates/MAP/other_rasters/accessibility/accessibility.5k.MEAN.tif'),
  Z('mastergrids/Other_Global_Covariates/Elevation/SRTM-Elevation/5km/Synoptic/SRTM_elevation.Synoptic.Overall.Data.5km.mean.tif'),
  Z('mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/5km/Synoptic/LST_Day_v6.Synoptic.Overall.SD.5km.mean.tif'),
  #Z('mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/TCB/5km/Synoptic/TCB.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/Other_Global_Covariates/NightTimeLights/VIIRS_DNB_Monthly/5km/Annual/VIIRS-SLC.2016.Annual.5km.MEDIAN.tif'),
  #Z('mastergrids/Other_Global_Covariates/UrbanAreas/Global_Urban_Footprint/From_86m/5km/Global_Urban_Footprint_5km_PropUrban.tif'),
  Z('mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/TCW_v6/5km/Synoptic/TCW_v6.Synoptic.Overall.mean.5km.mean.tif')
)


# load packages

## Spatial packages
library(raster)
library(maptools)
library(rgeos)

## dataframe packages
library(dplyr)
library(readr)
library(magrittr)
library(tidyr)

## For standardising prevalence
library(malariaAtlas)

## For inc prevalence conversions
library(GBDutils)

## Plotting packages
library(ggplot2)
library(cowplot)
theme_set(theme_minimal())

library(caret)



# load functions

source('collect_data.R')
source('process_data.R')
source('CombineRasters.R')
source('parallel-raster-extract.R')
source('build_inla_meshes.R')
source('fit_model.R')
source('run_cv.R')
source('random_crossvalidation_setup.R')
source('spatial_crossvalidation_setup.R')
source('plotting_functions.R')



set.seed(250918)

# load all data



# Read covariate rasters
covs_list <- lapply(cov_raster_paths, raster::raster)
crop_to <- find_smallest_extent(covs_list)
covs_cropped <- lapply(covs_list, function(x) crop(x, crop_to))
covs <- do.call(stack, CombineRasters(covs_cropped))






pr_pos_column = 'positive'
pr_n_column = 'examined'
pr_latlon = c('latitude', 'longitude')
pr_country = 'country'
pr_age_low = 'lower_age'
pr_age_high = 'upper_age'
standardisePR = c(2, 10)
roundPR = FALSE
standardisePars = 'Pf_Smith2007'
pr_min_year = 1990

pr <- readr::read_csv(PR_path, guess_max  = 1e5)

pr <- pr %>% filter(year_start >= pr_min_year)


pr_clean <- data_frame(
  prevalence = pull(pr, pr_pos_column) / pull(pr, pr_n_column),
  positive = pr %>% pull(pr_pos_column),
  examined = pr %>% pull(pr_n_column),
  latitude = pr %>% pull(pr_latlon[1]),
  longitude =  pr %>% pull(pr_latlon[2])
)
pr_clean <- pr_clean %>% mutate(prevalence = positive / examined)

if(!is.null(standardisePR)){
  prev_stand <- convertPrevalence(pr_clean$prevalence, 
                                  pr %>% pull(pr_age_low),
                                  pr %>% pull(pr_age_high),
                                  standardisePR[1],
                                  standardisePR[2],
                                  parameters = standardisePars
  )
  if(roundPR == TRUE){
    pr_clean$positive <- round(pr_clean$examined * prev_stand)
  } else {
    pr_clean$positive <- pr_clean$examined * prev_stand
  }
}  


coords <- pr_clean[, c(pr_latlon[2], pr_latlon[1])]

pr_extracted <- raster::extract(covs, SpatialPoints(coords))
names(pr_extracted) <- names(covs)


# ignoring GP_2013, remove NA rows.
missingpoints <- pr_extracted %>% complete.cases
pr_extracted <- pr_extracted[missingpoints, ]
pr_clean <- pr_clean[missingpoints, ]

#models <- fit_models()

#predictions <- predict_models()

m <- list()

y <- pr_clean$prevalence
partition <- createMultiFolds(y, k = 5, times = 1)

models <- c('enet', 'gbm', 'ranger', 'ppr', 'nnet')
tuneLength_vec <- c(10, 10, 10, 10, 10)
search_vec <- c('grid', 'random', 'random', 'grid', 'grid')

m[[1]] <- train(pr_extracted, y, 
                method = models[1],
                trControl = trainControl(index = partition, 
                                         returnData = TRUE,
                                         savePredictions = TRUE, 
                                         search = search_vec[1],
                                         predictionBounds = c(0, 1)),
                tuneLength = tuneLength_vec[1],
                weights = pr_clean$examined)




m[[2]] <- train(pr_extracted, y, 
                method = models[2],
                trControl = trainControl(index = partition, 
                                         returnData = TRUE,
                                         savePredictions = TRUE, 
                                         search = search_vec[2],
                                         predictionBounds = c(0, 1)),
                tuneLength = tuneLength_vec[2])






m[[3]] <- train(pr_extracted, y, 
                method = models[3],
                trControl = trainControl(index = partition, 
                                         returnData = TRUE,
                                         savePredictions = TRUE, 
                                         search = search_vec[3],
                                         predictionBounds = c(0, 1)),
                tuneLength = tuneLength_vec[3])



m[[4]] <- train(pr_extracted, y, 
                method = models[4],
                trControl = trainControl(index = partition, 
                                         returnData = TRUE,
                                         savePredictions = TRUE, 
                                         search = search_vec[4],
                                         predictionBounds = c(0, 1)),
                tuneLength = tuneLength_vec[4])


m[[5]] <- train(pr_extracted, y, 
                method = models[5],
                trControl = trainControl(index = partition, 
                                         returnData = TRUE,
                                         savePredictions = TRUE, 
                                         search = search_vec[5],
                                         predictionBounds = c(0, 1)),
                tuneLength = tuneLength_vec[5],
                linout = TRUE)






png('figs/enetopt_global.png')
print(plot(m[[1]]))
dev.off()

png('figs/ppropt_global.png')
print(plot(m[[4]]))
dev.off()

png('figs/nnetopt_global.png')
print(plot(m[[5]]))
dev.off()


p <- plotCV(m[[1]])
p + xlim(0, NA)
ggsave('figs/enet_obspred_global.png')


p <- plotCV(m[[2]])
p + xlim(0, NA)
ggsave('figs/gbm_obspred_global.png')

p <- plotCV(m[[3]])
p + xlim(0, NA)
ggsave('figs/ranger_obspred_global.png')

p <- plotCV(m[[4]])
p + xlim(0, NA)
ggsave('figs/ppr_obspred_global.png')



p <- plotCV(m[[5]])
p + xlim(0, NA)
ggsave('figs/nnet_obspred_global.png')


compare_models(m[[1]], m[[2]])
ggsave('figs/comp_enet_gbm_global.png')

compare_models(m[[3]], m[[2]])
ggsave('figs/comp_ranger_gbm_global.png')

compare_models(m[[4]], m[[2]])
ggsave('figs/comp_ppr_gbm_global.png')

compare_models(m[[5]], m[[2]])
ggsave('figs/comp_nnet_gbm_global.png')

save(m, file = 'model_outputs/global_ml_global.RData')


extent_col <- c(-85, -60, -10, 20)
covs_crop_col <- crop(covs, extent_col)
covs_col_mat <- getValues(covs_crop_col)

pred <- matrix(NA, nrow = nrow(covs_col_mat), ncol = length(m))

nas <- complete.cases(covs_col_mat)


pred[nas, ] <- predict(m, newdata = covs_col_mat[nas, ], na.action = na.pass) %>% do.call(cbind, .)

r.pts <- rasterToPoints(covs_crop_col, spatial = TRUE)


pred_rast_col <- rasterFromXYZ(cbind(r.pts@coords, pred))
pred_rast_col[pred_rast_col < 0] <- 0
pred_rast_col_inc <- calc(pred_rast_col, PrevIncConversion)
names(pred_rast_col_inc) <- sapply(m, function(x) x$method)########
# master script for point Vs polygon Vs joint analysis
# Tim Lucas
# 2018-05-30
###########

if(Sys.info()["user"] != 'anita'){
  setwd('~/timz/timothy/polygon_ml')
} else {
  #setwd('~/Z/users/anita/point_polygon_join_comparison_analysis')
}

source("setUserInfo.R")

# define paths

PR_path <- Z('GBD2017/Processing/Stages/04b_PR_DB_Import_Export/Verified_Outputs/2018_02_15/pfpr_dump.csv')
API_path <- Z('GBD2017/Processing/Stages/04c_API_Data_Export/Checkpoint_Outputs/subnational.csv')
pop_path <- Z('GBD2017/Processing/Stages/03_Muster_Population_Figures/Verified_Outputs/Output_Pop_At_Risk_Pf_5K/ihme_corrected_frankenpop_All_Ages_3_2015_at_risk_pf.tif')
shapefile_path <- Z('master_geometries/Admin_Units/Global/GBD/GBD2017_MAP/GBD2017_MAP_MG_5K/')


cov_raster_paths <- c(
  Z('mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/5km/Synoptic/LST_Day_v6.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/EVI_v6/5km/Synoptic/EVI_v6.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Synoptic/TSI-Martens2-Pf.Synoptic.Overall.Mean.5km.Data.tif'),
  Z('GBD2017/Processing/Static_Covariates/MAP/other_rasters/accessibility/accessibility.5k.MEAN.tif'),
  Z('mastergrids/Other_Global_Covariates/Elevation/SRTM-Elevation/5km/Synoptic/SRTM_elevation.Synoptic.Overall.Data.5km.mean.tif'),
  Z('mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/5km/Synoptic/LST_Day_v6.Synoptic.Overall.SD.5km.mean.tif'),
  #Z('mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/TCB/5km/Synoptic/TCB.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/Other_Global_Covariates/NightTimeLights/VIIRS_DNB_Monthly/5km/Annual/VIIRS-SLC.2016.Annual.5km.MEDIAN.tif'),
  #Z('mastergrids/Other_Global_Covariates/UrbanAreas/Global_Urban_Footprint/From_86m/5km/Global_Urban_Footprint_5km_PropUrban.tif'),
  Z('mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/TCW_v6/5km/Synoptic/TCW_v6.Synoptic.Overall.mean.5km.mean.tif')
)


# load packages

## Spatial packages
library(raster)
library(maptools)
library(rgeos)

## dataframe packages
library(dplyr)
library(readr)
library(magrittr)
library(tidyr)

## For standardising prevalence
library(malariaAtlas)

## For inc prevalence conversions
library(GBDutils)

## Plotting packages
library(ggplot2)
library(cowplot)
theme_set(theme_minimal())

library(caret)



# load functions

source('collect_data.R')
source('process_data.R')
source('CombineRasters.R')
source('parallel-raster-extract.R')
source('build_inla_meshes.R')
source('fit_model.R')
source('run_cv.R')
source('random_crossvalidation_setup.R')
source('spatial_crossvalidation_setup.R')
source('plotting_functions.R')



set.seed(250918)

# load all data



# Read covariate rasters
covs_list <- lapply(cov_raster_paths, raster::raster)
crop_to <- find_smallest_extent(covs_list)
covs_cropped <- lapply(covs_list, function(x) crop(x, crop_to))
covs <- do.call(stack, CombineRasters(covs_cropped))






pr_pos_column = 'positive'
pr_n_column = 'examined'
pr_latlon = c('latitude', 'longitude')
pr_country = 'country'
pr_age_low = 'lower_age'
pr_age_high = 'upper_age'
standardisePR = c(2, 10)
roundPR = FALSE
standardisePars = 'Pf_Smith2007'
pr_min_year = 1990

pr <- readr::read_csv(PR_path, guess_max  = 1e5)

pr_region <- 'SouthAmerica'
if(pr_region == 'country'){
  usecountries <- find_country_from_iso3(useiso3, api_full$iso3, api_full$country_name)
} else if(pr_region == 'SouthAmerica'){
  usecountries <- c("Colombia", "Brazil", "Venezuela", "Peru", "Suriname", "Bolivia")
} else if(pr_region == 'all'){
  usecountries <- unique(pr$country)
}
pr <- pr %>% filter(country %in% usecountries, year_start >= pr_min_year)


pr_clean <- data_frame(
  prevalence = pull(pr, pr_pos_column) / pull(pr, pr_n_column),
  positive = pr %>% pull(pr_pos_column),
  examined = pr %>% pull(pr_n_column),
  latitude = pr %>% pull(pr_latlon[1]),
  longitude =  pr %>% pull(pr_latlon[2])
)
pr_clean <- pr_clean %>% mutate(prevalence = positive / examined)

if(!is.null(standardisePR)){
  prev_stand <- convertPrevalence(pr_clean$prevalence, 
                                  pr %>% pull(pr_age_low),
                                  pr %>% pull(pr_age_high),
                                  standardisePR[1],
                                  standardisePR[2],
                                  parameters = standardisePars
  )
  if(roundPR == TRUE){
    pr_clean$positive <- round(pr_clean$examined * prev_stand)
  } else {
    pr_clean$positive <- pr_clean$examined * prev_stand
  }
}  


coords <- pr_clean[, c(pr_latlon[2], pr_latlon[1])]

pr_extracted <- raster::extract(covs, SpatialPoints(coords))
names(pr_extracted) <- names(covs)


# ignoring GP_2013, remove NA rows.
missingpoints <- pr_extracted %>% complete.cases
pr_extracted <- pr_extracted[missingpoints, ]
pr_clean <- pr_clean[missingpoints, ]

#models <- fit_models()

#predictions <- predict_models()

m <- list()

y <- pr_clean$prevalence
partition <- createMultiFolds(y, k = 5, times = 1)

models <- c('enet', 'gbm', 'ranger', 'ppr', 'nnet')
tuneLength_vec <- c(10, 10, 10, 10, 10)
search_vec <- c('grid', 'random', 'random', 'grid', 'grid')

m[[1]] <- train(pr_extracted, y, 
                method = models[1],
                trControl = trainControl(index = partition, 
                                         returnData = TRUE,
                                         savePredictions = TRUE, 
                                         search = search_vec[1],
                                         predictionBounds = c(0, 1)),
                tuneLength = tuneLength_vec[1],
                weights = pr_clean$examined)




m[[2]] <- train(pr_extracted, y, 
                method = models[2],
                trControl = trainControl(index = partition, 
                                         returnData = TRUE,
                                         savePredictions = TRUE, 
                                         search = search_vec[2],
                                         predictionBounds = c(0, 1)),
                tuneLength = tuneLength_vec[2])






m[[3]] <- train(pr_extracted, y, 
                method = models[3],
                trControl = trainControl(index = partition, 
                                         returnData = TRUE,
                                         savePredictions = TRUE, 
                                         search = search_vec[3],
                                         predictionBounds = c(0, 1)),
                tuneLength = tuneLength_vec[3])



m[[4]] <- train(pr_extracted, y, 
                method = models[4],
                trControl = trainControl(index = partition, 
                                         returnData = TRUE,
                                         savePredictions = TRUE, 
                                         search = search_vec[4],
                                         predictionBounds = c(0, 1)),
                tuneLength = tuneLength_vec[4])


m[[5]] <- train(pr_extracted, y, 
                method = models[5],
                trControl = trainControl(index = partition, 
                                         returnData = TRUE,
                                         savePredictions = TRUE, 
                                         search = search_vec[5],
                                         predictionBounds = c(0, 1)),
                tuneLength = tuneLength_vec[5],
                linout = TRUE)






png('figs/enetopt_global.png')
print(plot(m[[1]]))
dev.off()

png('figs/ppropt_global.png')
print(plot(m[[4]]))
dev.off()

png('figs/nnetopt_global.png')
print(plot(m[[5]]))
dev.off()


p <- plotCV(m[[1]])
p + xlim(0, NA)
ggsave('figs/enet_obspred_global.png')


p <- plotCV(m[[2]])
p + xlim(0, NA)
ggsave('figs/gbm_obspred_global.png')

p <- plotCV(m[[3]])
p + xlim(0, NA)
ggsave('figs/ranger_obspred_global.png')

p <- plotCV(m[[4]])
p + xlim(0, NA)
ggsave('figs/ppr_obspred_global.png')



p <- plotCV(m[[5]])
p + xlim(0, NA)
ggsave('figs/nnet_obspred_global.png')


compare_models(m[[1]], m[[2]])
ggsave('figs/comp_enet_gbm_global.png')

compare_models(m[[3]], m[[2]])
ggsave('figs/comp_ranger_gbm_global.png')

compare_models(m[[4]], m[[2]])
ggsave('figs/comp_ppr_gbm_global.png')

compare_models(m[[5]], m[[2]])
ggsave('figs/comp_nnet_gbm_global.png')

save(m, file = 'model_outputs/global_ml_global.RData')


extent_col <- c(-85, -60, -10, 20)
covs_crop_col <- crop(covs, extent_col)
covs_col_mat <- getValues(covs_crop_col)

pred <- matrix(NA, nrow = nrow(covs_col_mat), ncol = length(m))

nas <- complete.cases(covs_col_mat)


pred[nas, ] <- predict(m, newdata = covs_col_mat[nas, ], na.action = na.pass) %>% do.call(cbind, .)

r.pts <- rasterToPoints(covs_crop_col, spatial = TRUE)

########
# master script for point Vs polygon Vs joint analysis
# Tim Lucas
# 2018-05-30
###########

if(Sys.info()["user"] != 'anita'){
  setwd('~/timz/timothy/polygon_ml')
} else {
  #setwd('~/Z/users/anita/point_polygon_join_comparison_analysis')
}

source("setUserInfo.R")

# define paths

PR_path <- Z('GBD2017/Processing/Stages/04b_PR_DB_Import_Export/Verified_Outputs/2018_02_15/pfpr_dump.csv')
API_path <- Z('GBD2017/Processing/Stages/04c_API_Data_Export/Checkpoint_Outputs/subnational.csv')
pop_path <- Z('GBD2017/Processing/Stages/03_Muster_Population_Figures/Verified_Outputs/Output_Pop_At_Risk_Pf_5K/ihme_corrected_frankenpop_All_Ages_3_2015_at_risk_pf.tif')
shapefile_path <- Z('master_geometries/Admin_Units/Global/GBD/GBD2017_MAP/GBD2017_MAP_MG_5K/')


cov_raster_paths <- c(
  Z('mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/5km/Synoptic/LST_Day_v6.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/EVI_v6/5km/Synoptic/EVI_v6.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Synoptic/TSI-Martens2-Pf.Synoptic.Overall.Mean.5km.Data.tif'),
  Z('GBD2017/Processing/Static_Covariates/MAP/other_rasters/accessibility/accessibility.5k.MEAN.tif'),
  Z('mastergrids/Other_Global_Covariates/Elevation/SRTM-Elevation/5km/Synoptic/SRTM_elevation.Synoptic.Overall.Data.5km.mean.tif'),
  Z('mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/5km/Synoptic/LST_Day_v6.Synoptic.Overall.SD.5km.mean.tif'),
  #Z('mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/TCB/5km/Synoptic/TCB.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/Other_Global_Covariates/NightTimeLights/VIIRS_DNB_Monthly/5km/Annual/VIIRS-SLC.2016.Annual.5km.MEDIAN.tif'),
  #Z('mastergrids/Other_Global_Covariates/UrbanAreas/Global_Urban_Footprint/From_86m/5km/Global_Urban_Footprint_5km_PropUrban.tif'),
  Z('mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/TCW_v6/5km/Synoptic/TCW_v6.Synoptic.Overall.mean.5km.mean.tif')
)


# load packages

## Spatial packages
library(raster)
library(maptools)
library(rgeos)

## dataframe packages
library(dplyr)
library(readr)
library(magrittr)
library(tidyr)

## For standardising prevalence
library(malariaAtlas)

## For inc prevalence conversions
library(GBDutils)

## Plotting packages
library(ggplot2)
library(cowplot)
theme_set(theme_minimal())

library(caret)



# load functions

source('collect_data.R')
source('process_data.R')
source('CombineRasters.R')
source('parallel-raster-extract.R')
source('build_inla_meshes.R')
source('fit_model.R')
source('run_cv.R')
source('random_crossvalidation_setup.R')
source('spatial_crossvalidation_setup.R')
source('plotting_functions.R')



set.seed(250918)

# load all data



# Read covariate rasters
covs_list <- lapply(cov_raster_paths, raster::raster)
crop_to <- find_smallest_extent(covs_list)
covs_cropped <- lapply(covs_list, function(x) crop(x, crop_to))
covs <- do.call(stack, CombineRasters(covs_cropped))






pr_pos_column = 'positive'
pr_n_column = 'examined'
pr_latlon = c('latitude', 'longitude')
pr_country = 'country'
pr_age_low = 'lower_age'
pr_age_high = 'upper_age'
standardisePR = c(2, 10)
roundPR = FALSE
standardisePars = 'Pf_Smith2007'
pr_min_year = 1990

pr <- readr::read_csv(PR_path, guess_max  = 1e5)

pr_region <- 'SouthAmerica'
if(pr_region == 'country'){
  usecountries <- find_country_from_iso3(useiso3, api_full$iso3, api_full$country_name)
} else if(pr_region == 'SouthAmerica'){
  usecountries <- c("Colombia", "Brazil", "Venezuela", "Peru", "Suriname", "Bolivia")
} else if(pr_region == 'all'){
  usecountries <- unique(pr$country)
}
pr <- pr %>% filter(country %in% usecountries, year_start >= pr_min_year)


pr_clean <- data_frame(
  prevalence = pull(pr, pr_pos_column) / pull(pr, pr_n_column),
  positive = pr %>% pull(pr_pos_column),
  examined = pr %>% pull(pr_n_column),
  latitude = pr %>% pull(pr_latlon[1]),
  longitude =  pr %>% pull(pr_latlon[2])
)
pr_clean <- pr_clean %>% mutate(prevalence = positive / examined)

if(!is.null(standardisePR)){
  prev_stand <- convertPrevalence(pr_clean$prevalence, 
                                  pr %>% pull(pr_age_low),
                                  pr %>% pull(pr_age_high),
                                  standardisePR[1],
                                  standardisePR[2],
                                  parameters = standardisePars
  )
  if(roundPR == TRUE){
    pr_clean$positive <- round(pr_clean$examined * prev_stand)
  } else {
    pr_clean$positive <- pr_clean$examined * prev_stand
  }
}  


coords <- pr_clean[, c(pr_latlon[2], pr_latlon[1])]

pr_extracted <- raster::extract(covs, SpatialPoints(coords))
names(pr_extracted) <- names(covs)


# ignoring GP_2013, remove NA rows.
missingpoints <- pr_extracted %>% complete.cases
pr_extracted <- pr_extracted[missingpoints, ]
pr_clean <- pr_clean[missingpoints, ]

#models <- fit_models()

#predictions <- predict_models()

m <- list()

y <- pr_clean$prevalence
partition <- createMultiFolds(y, k = 3, times = 1)

models <- c('enet', 'xgbTree', 'ranger', 'ppr', 'nnet')
tuneLength_vec <- c(10, 2, 2, 10, 2)
search_vec <- c('grid', 'random', 'random', 'grid', 'grid')

m[[1]] <- train(pr_extracted, y, 
                method = models[1],
                trControl = trainControl(index = partition, 
                                         returnData = TRUE,
                                         savePredictions = TRUE, 
                                         search = search_vec[1],
                                         predictionBounds = c(0, 1)),
                tuneLength = tuneLength_vec[1],
                weights = pr_clean$examined)




m[[2]] <- train(pr_extracted, y, 
                method = models[2],
                trControl = trainControl(index = partition, 
                                         returnData = TRUE,
                                         savePredictions = TRUE, 
                                         search = search_vec[2],
                                         predictionBounds = c(0, 1)),
                tuneLength = tuneLength_vec[2],
                base_score = mean(y))






m[[3]] <- train(pr_extracted, y, 
                method = models[3],
                trControl = trainControl(index = partition, 
                                         returnData = TRUE,
                                         savePredictions = TRUE, 
                                         search = search_vec[3],
                                         predictionBounds = c(0, 1)),
                tuneLength = tuneLength_vec[3])



m[[4]] <- train(pr_extracted, y, 
                method = models[4],
                trControl = trainControl(index = partition, 
                                         returnData = TRUE,
                                         savePredictions = TRUE, 
                                         search = search_vec[4],
                                         predictionBounds = c(0, 1)),
                tuneLength = tuneLength_vec[4])


m[[5]] <- train(pr_extracted, y, 
                method = models[5],
                trControl = trainControl(index = partition, 
                                         returnData = TRUE,
                                         savePredictions = TRUE, 
                                         search = search_vec[5],
                                         predictionBounds = c(0, 1)),
                tuneLength = tuneLength_vec[5],
                linout = TRUE)






png('figs/enetopt_global.png')
print(plot(m[[1]]))
dev.off()

png('figs/ppropt_global.png')
print(plot(m[[4]]))
dev.off()

png('figs/nnetopt_global.png')
print(plot(m[[5]]))
dev.off()


p <- plotCV(m[[1]])
p + xlim(0, NA)
ggsave('figs/enet_obspred_global.png')


p <- plotCV(m[[2]])
p + xlim(0, NA)
ggsave('figs/gbm_obspred_global.png')

p <- plotCV(m[[3]])
p + xlim(0, NA)
ggsave('figs/ranger_obspred_global.png')

p <- plotCV(m[[4]])
p + xlim(0, NA)
ggsave('figs/ppr_obspred_global.png')



p <- plotCV(m[[5]])
p + xlim(0, NA)
ggsave('figs/nnet_obspred_global.png')


compare_models(m[[1]], m[[2]])
ggsave('figs/comp_enet_gbm_global.png')

compare_models(m[[3]], m[[2]])
ggsave('figs/comp_ranger_gbm_global.png')

compare_models(m[[4]], m[[2]])
ggsave('figs/comp_ppr_gbm_global.png')

compare_models(m[[5]], m[[2]])
ggsave('figs/comp_nnet_gbm_global.png')


compare_models(m[[1]], m[[4]])
ggsave('figs/comp_enet_ppr_global.png')

save(m, file = 'model_outputs/global_ml_global.RData')







