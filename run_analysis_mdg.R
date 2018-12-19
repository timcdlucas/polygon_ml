########
# master script for ml on points vs polygon only analysis
# Tim Lucas
# 2018-08-20
###########

if(Sys.info()["user"] != 'anita'){
  setwd('~/timz/timothy/polygon_ml_wsc')
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

ml_local_raster_paths <- c(
  'model_outputs/ml_pred_rasters/madagascar_mdg_enet.tif',
  'model_outputs/ml_pred_rasters/madagascar_mdg_gbm.tif',
  'model_outputs/ml_pred_rasters/madagascar_mdg_nnet.tif',
  'model_outputs/ml_pred_rasters/madagascar_mdg_ppr.tif',
  'model_outputs/ml_pred_rasters/madagascar_mdg_ranger.tif'
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

#library(caret)
##  Modelling packages
library(TMB)
#library(stantmb)
library(INLA)
if(Sys.info()["sysname"] != 'Windows'){
  message('using INLA unix workaround')
  INLA:::inla.dynload.workaround()
} else {
  message('Not using INLA unix workaround. Expect you are using winows.')
}
library(INLAutils)
library(sparseMVN)

# Parallel processing
library(foreach)
library(doParallel)


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

# Compile the model
compile("joint_model.cpp")

set.seed(180530)

# load all data
# Perhaps read all raster years and have a step in process data to choose the right one? Or something. Kinda annoying.
data <- load_data(PR_path, 
                  API_path, 
                  pop_path, 
                  cov_raster_paths, 
                  shapefile_path, 
                  shapefile_pattern = '.shp$', 
                  useiso3 = 'MDG', 
                  admin_unit_level = 'ADMIN3',
                  pr_country = 'country',
                  api_year = 2013)

data_ml_cov <- load_data(PR_path, 
                         API_path, 
                         pop_path, 
                         ml_local_raster_paths, 
                         shapefile_path, 
                         shapefile_pattern = '.shp$', 
                         useiso3 = 'MDG', 
                         admin_unit_level = 'ADMIN3',
                         pr_country = 'country',
                         api_year = 2013)

data_all_cov <- load_data(PR_path, 
                          API_path, 
                          pop_path, 
                          c(cov_raster_paths, ml_local_raster_paths), 
                          shapefile_path, 
                          shapefile_pattern = '.shp$', 
                          useiso3 = 'MDG', 
                          admin_unit_level = 'ADMIN3',
                          pr_country = 'country',
                          api_year = 2013)


# pre analysis

data_mdg_cov <- process_data(
  binomial_positive = data$pr$positive,
  binomial_n = data$pr$examined,
  coords = data$pr[, c('longitude', 'latitude')],
  polygon_response = data$api$api_mean,
  polygon_population = data$api$population,
  shapefile_id = data$api$shapefile_id,
  shps_id_column = 'area_id',
  shapefiles = data$shapefiles,
  pop_raster = data$pop,
  cov_rasters = data$covs,
  useiso3 = 'MDG',
  transform = c(4:7))

save(data_mdg_cov, file = 'model_outputs/mdg_cov_data.RData')

data_mdg_ml <- process_data(
  binomial_positive = data_ml_cov$pr$positive,
  binomial_n = data_ml_cov$pr$examined,
  coords = data_ml_cov$pr[, c('longitude', 'latitude')],
  polygon_response = data_ml_cov$api$api_mean,
  polygon_population = data_ml_cov$api$population,
  shapefile_id = data_ml_cov$api$shapefile_id,
  shps_id_column = 'area_id',
  shapefiles = data_ml_cov$shapefiles,
  pop_raster = data_ml_cov$pop,
  cov_rasters = data_ml_cov$covs,
  useiso3 = 'MDG',
  transform = NULL)
save(data_mdg_ml, file = 'model_outputs/mdg_ml_data.RData')

data_mdg_all <- process_data(
  binomial_positive = data_all_cov$pr$positive,
  binomial_n = data_all_cov$pr$examined,
  coords = data_all_cov$pr[, c('longitude', 'latitude')],
  polygon_response = data_all_cov$api$api_mean,
  polygon_population = data_all_cov$api$population,
  shapefile_id = data_all_cov$api$shapefile_id,
  shps_id_column = 'area_id',
  shapefiles = data_all_cov$shapefiles,
  pop_raster = data_all_cov$pop,
  cov_rasters = data_all_cov$covs,
  useiso3 = 'MDG',
  transform = c(4:7))

save(data_mdg_all, file = 'model_outputs/mdg_all_data.RData')


autoplot(data_mdg_cov, pr_limits = c(0, 0.3))
autoplot(data_mdg_cov, pr_limits = c(0, 0.3), trans = 'log1p')

ggsave('figs/mdg_input_data.png')



mesh_mdg <- build_mesh(data_mdg_cov, mesh.args = list(max.edge = c(0.4, 4), cut = 0.4))
autoplot(mesh_mdg)
save(mesh_mdg, file = 'model_outputs/mdg_mesh.RData')



# Define cross validation strategies
data_cv1_mdg <- cv_random_folds(data_mdg_cov, k = 6)
data_cv1_mdg_ml <- cv_random_folds(data_mdg_ml, k = 6, 
                                   polygon_folds = attr(data_cv1_mdg, 'polygon_folds'),
                                   pr_folds = attr(data_cv1_mdg, 'pr_folds'))
data_cv1_mdg_all <- cv_random_folds(data_mdg_all, k = 6, 
                                    polygon_folds = attr(data_cv1_mdg, 'polygon_folds'),
                                    pr_folds = attr(data_cv1_mdg, 'pr_folds'))


autoplot(data_cv1_mdg, jitter = 0)
autoplot(data_cv1_mdg_ml, jitter = 0)

ggsave('figs/mdg_cv_random.png')
save(data_cv1_mdg, file = 'model_outputs/mdg_cv_1.RData')


# Spatial
data_cv2_mdg <- cv_spatial_folds(data_mdg_cov, k = 3)
data_cv2_mdg_ml <- cv_spatial_folds(data_mdg_ml, k = 3, 
                                    polygon_folds = attr(data_cv2_mdg, 'polygon_folds'),
                                    pr_folds = attr(data_cv2_mdg, 'pr_folds'))
data_cv2_mdg_all <- cv_spatial_folds(data_mdg_all, k = 3, 
                                     polygon_folds = attr(data_cv2_mdg, 'polygon_folds'),
                                     pr_folds = attr(data_cv2_mdg, 'pr_folds'))

autoplot(data_cv2_mdg, jitter = 0)
autoplot(data_cv2_mdg_ml, jitter = 0)

ggsave('figs/mdg_cv_spatial2.png')
save(data_cv2_mdg, file = 'model_outputs/mdg_cv_2.RData')


#autoplot(data_cv1_mdg[[1]]$train, pr_limits = c(0, 0.3))

use_points <- 0
use_polygons <- 1
# run models
# Run full model to get a handle on things.

arg_list <- list(prior_rho_min = 1, # 
                 prior_rho_prob = 0.00001, # Want p(rho < 3) = 0.0001
                 prior_sigma_max = 1, # Want p(sd > 1) = 0.0001 (would explain most of prev). 
                 prior_sigma_prob = 0.00001,
                 prior_iideffect_sd_max = 0.05, 
                 # The difference between m_low_pf and LCI(pois(m_mean_pf)), then converted to inc rate, then to prev ranges around 0-0.025. 
                 # The 0.975 quantile of that (two sided) distribution is 0.005 prevalence. 
                 # To explain 0.005 prevalence, we need a norm of 0.05. Fin.
                 prior_iideffect_sd_prob = 0.0000001, # Made this stronger because too much iid.
                 prior_iideffect_pr_sd_max = 0.3, # Max difference between PR points within a cell (with n > 500)
                 prior_iideffect_pr_sd_prob = 0.0000001,
                 priormean_intercept = -2,
                 priorsd_intercept = 2,  # Indonesia has prev lowish. But want intercept to take whatever value it likes.
                 priormean_slope = 0, 
                 priorsd_slope = 0.4, # Explains between 0.004 and 0.27 prevalence. 1 covariate shouldn't explain between 0 and 0.6 (range of prev).
                 use_polygons = use_polygons,
                 use_points = use_points)

if(FALSE){
  full_model <- fit_model(data_mdg_cov, mesh_mdg, its = 1000, model.args = arg_list)
  autoplot(full_model)
  
  png('figs/full_model_covs_in_sample_map.png')
  plot(full_model, layer = 'api')
  dev.off()
  
  png('figs/full_model_covs_covariates_in_sample_map.png')
  plot(full_model, layer = 'covariates')
  dev.off()
  
  png('figs/full_model_covs_field_in_sample_map.png')
  plot(full_model, layer = 'field')
  dev.off()
  
  png('figs/full_model_covs_in_sample_map_log.png')
  full_model$predictions$api %>% log10 %>% plot
  dev.off()
  
  
  
  in_sample <- cv_performance(predictions = full_model$predictions, 
                              holdout = data_mdg_cov,
                              model_params = full_model$model, 
                              CI = 0.8,
                              use_points = use_points)
  autoplot(in_sample, CI = TRUE)
  autoplot(in_sample, trans = 'log1p', CI = TRUE)
  ggsave('figs/mdg_full_model_covs_in_sample.png')
  
  save(full_model, file = 'model_outputs/full_model_covs_mdg.RData')
  
  
  
  
  full_model_ml <- fit_model(data_ml_cov, mesh_mdg, its = 1000, model.args = arg_list)
  autoplot(full_model_ml)
  
  png('figs/full_model_ml_in_sample_map.png')
  plot(full_model_ml, layer = 'api')
  dev.off()
  
  png('figs/full_model_ml_covariates_in_sample_map.png')
  plot(full_model_ml, layer = 'covariates')
  dev.off()
  
  png('figs/full_model_ml_field_in_sample_map.png')
  plot(full_model_ml, layer = 'field')
  dev.off()
  
  png('figs/full_model_ml_in_sample_map_log.png')
  full_model_ml$predictions$api %>% log10 %>% plot
  dev.off()
  
  
  
  in_sample_ml <- cv_performance(predictions = full_model_ml$predictions, 
                                 holdout = data_ml_cov,
                                 model_params = full_model_ml$model, 
                                 CI = 0.8,
                                 use_points = use_points)
  autoplot(in_sample_ml, CI = TRUE)
  autoplot(in_sample_ml, trans = 'log1p', CI = TRUE)
  ggsave('figs/mdg_full_model_ml_in_sample.png')
  
  save(full_model_ml, file = 'model_outputs/full_model_ml_mdg.RData')
  
  
  
  
  
  full_model_all <- fit_model(data_mdg_all, mesh_mdg, its = 1000, model.args = arg_list)
  autoplot(full_model_all)
  
  png('figs/full_model_all_in_sample_map.png')
  plot(full_model_all, layer = 'api')
  dev.off()
  
  png('figs/full_model_all_in_sample_map_log.png')
  full_model_all$predictions$api %>% log10 %>% plot
  dev.off()
  
  
  png('figs/full_model_all_covariates_in_sample_map.png')
  plot(full_model_all, layer = 'covariates')
  dev.off()
  
  png('figs/full_mode_alll_field_in_sample_map.png')
  plot(full_model_all, layer = 'field')
  dev.off()
  
  in_sample_all <- cv_performance(predictions = full_model_all$predictions, 
                                  holdout = data_mdg_all,
                                  model_params = full_model_all$model, 
                                  CI = 0.8,
                                  use_points = use_points)
  autoplot(in_sample_all, CI = TRUE)
  autoplot(in_sample_all, trans = 'log1p', CI = TRUE)
  ggsave('figs/mdg_full_model_all_in_sample.png')
  
  save(in_sample_all, file = 'model_outputs/full_model_all_mdg_all.RData')
  
  
}


# Run 3 x models on cv1.
cat('Start cv1 model 1\n')

cv1_output1 <- run_cv(data_cv1_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, parallel_delay = 20, cores = 3)
obspred_map(data_cv1_mdg, cv1_output1, column = FALSE, mask = TRUE)
ggsave('figs/mdg_covs_only_obspred_map.png')
obspred_map(data_cv1_mdg, cv1_output1, trans = 'log10', column = FALSE, mask = TRUE)
ggsave('figs/mdg_covs_oonly_obspred_map_log.png')
autoplot(cv1_output1, type = 'obs_preds', CI = F)
ggsave('figs/mdg_covs_oonly_obspred.png')
autoplot(cv1_output1, type = 'obs_preds', CI = FALSE, tran = 'log1p')
ggsave('figs/mdg_covs_only_obspred_log.png')

cat('Start cv1 model 2\n')

cv1_output2 <- run_cv(data_cv1_mdg_ml, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, parallel_delay = 40, cores = 3)
obspred_map(data_cv1_mdg, cv1_output2, column = FALSE)
ggsave('figs/mdg_ml_only_obspred_map.png')
obspred_map(data_cv1_mdg, cv1_output2, trans = 'log10', column = FALSE)
ggsave('figs/mdg_ml_only_obspred_map_log.png')
autoplot(cv1_output2, type = 'obs_preds', CI = FALSE)
ggsave('figs/mdg_ml_only_obspred.png')
autoplot(cv1_output2, type = 'obs_preds', CI = FALSE, tran = 'log1p')
ggsave('figs/mdg_ml_only_obspred_log.png')

cat('Start cv1 model 3\n')

cv1_output3 <- run_cv(data_cv1_mdg_all, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, parallel_delay = 20, cores = 3)
obspred_map(data_cv1_mdg, cv1_output3, column = FALSE)
ggsave('figs/mdg_all_obspred_map.png')
obspred_map(data_cv1_mdg, cv1_output3, trans = 'log10', column = FALSE)
ggsave('figs/mdg_all_obspred_map_log.png')
autoplot(cv1_output3, type = 'obs_preds', CI = FALSE)
ggsave('figs/mdg_all_obspred.png')
autoplot(cv1_output3, type = 'obs_preds', CI = FALSE, tran = 'log1p')
ggsave('figs/mdg_all_obspred_log.png')


save(cv1_output1, file = 'model_outputs/mdg_covs_cv_1.RData')
save(cv1_output2, file = 'model_outputs/mdg_ml_cv_1.RData')
save(cv1_output3, file = 'model_outputs/mdg_all_cv_1.RData')

cv1_output1$summary$polygon_metrics
cv1_output2$summary$polygon_metrics
cv1_output3$summary$polygon_metrics

cv1_output1$summary$pr_metrics
cv1_output2$summary$pr_metrics
cv1_output3$summary$pr_metrics




# Run 3 x models on cv2.
cat('Start cv2 model 1')

cv2_output1 <- run_cv(data_cv2_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, parallel_delay = 40, cores = 3)
obspred_map(data_cv2_mdg, cv2_output1, column = FALSE, mask = TRUE)
ggsave('figs/mdg_covs_only_obspred_map2.png')
obspred_map(data_cv2_mdg, cv2_output1, trans = 'log10', column = FALSE, mask = TRUE)
ggsave('figs/mdg_covs_only_obspred_map_log2.png')
autoplot(cv2_output1, type = 'obs_preds', CI = FALSE)
ggsave('figs/mdg_covs_only_obspred2.png')
autoplot(cv2_output1, type = 'obs_preds', CI = FALSE, tran = 'log1p')
ggsave('figs/mdg_covs_only_obspred_log2.png')


cat('Start cv2 model 2')

cv2_output2 <- run_cv(data_cv2_mdg_ml, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, parallel_delay = 30, cores = 3)
obspred_map(data_cv2_mdg, cv2_output2, column = FALSE, mask = TRUE)
ggsave('figs/mdg_ml_only_obspred_map2.png')
obspred_map(data_cv2_mdg, cv2_output2, trans = 'log10', column = FALSE, mask = TRUE)
ggsave('figs/mdg_ml_only_obspred_map_log2.png')
autoplot(cv2_output2, type = 'obs_preds', CI = FALSE)
ggsave('figs/mdg_ml_only_obspred2.png')
autoplot(cv2_output2, type = 'obs_preds', CI = FALSE, tran = 'log1p')
ggsave('figs/mdg_ml_only_obspred_log2.png')


cat('Start cv2 model 3')

cv2_output3 <- run_cv(data_cv2_mdg_all, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, parallel_delay = 0, cores = 3)
obspred_map(data_cv2_mdg, cv2_output3, column = FALSE, mask = TRUE)
ggsave('figs/mdg_all_obspred_map2.png')
obspred_map(data_cv2_mdg, cv2_output3, trans = 'log10', column = FALSE, mask = TRUE)
ggsave('figs/mdg_all_obspred_map_log2.png')
autoplot(cv2_output3, type = 'obs_preds', CI = FALSE)
ggsave('figs/mdg_all_obspred2.png')
autoplot(cv2_output3, type = 'obs_preds', CI = FALSE, tran = 'log1p')
ggsave('figs/mdg_all_only_obspred_log2.png')



save(cv2_output1, file = 'model_outputs/mdg_covs_cv_2.RData')
save(cv2_output2, file = 'model_outputs/mdg_ml_cv_2.RData')
save(cv2_output3, file = 'model_outputs/mdg_all_cv_2.RData')

cv2_output1$summary$polygon_metrics
cv2_output2$summary$polygon_metrics
cv2_output3$summary$polygon_metrics

cv2_output1$summary$pr_metrics
cv2_output2$summary$pr_metrics
cv2_output3$summary$pr_metrics







