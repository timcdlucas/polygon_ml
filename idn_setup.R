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
  Z('mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/5km/Synoptic/LST_Day.Synoptic.Overall.mean.5km.mean.tif'),
  #Z('mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/EVI/5km/Synoptic/EVI.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Synoptic/TSI-Martens2-Pf.Synoptic.Overall.Mean.5km.Data.tif'),
  Z('GBD2017/Processing/Static_Covariates/MAP/other_rasters/accessibility/accessibility.5k.MEAN.tif'),
  Z('mastergrids/Other_Global_Covariates/Elevation/SRTM-Elevation/5km/Synoptic/SRTM_elevation.Synoptic.Overall.Data.5km.mean.tif'),
  Z('mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/5km/Synoptic/LST_Day.Synoptic.Overall.SD.5km.mean.tif'),
  #Z('mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/TCB/5km/Synoptic/TCB.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/Other_Global_Covariates/NightTimeLights/VIIRS_DNB_Monthly/5km/Annual/VIIRS-SLC.2016.Annual.5km.MEDIAN.tif'),
  #Z('mastergrids/Other_Global_Covariates/UrbanAreas/Global_Urban_Footprint/From_86m/5km/Global_Urban_Footprint_5km_PropUrban.tif'),
  Z('mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/TCW/5km/Synoptic/TCW.Synoptic.Overall.mean.5km.mean.tif')
)

ml_local_raster_paths <- c(
  'model_outputs/ml_pred_rasters/south_asia_idn_enet.tif',
  'model_outputs/ml_pred_rasters/south_asia_idn_ppr.tif',
  'model_outputs/ml_pred_rasters/south_asia_idn_xgbTree.tif',
  'model_outputs/ml_pred_rasters/south_asia_idn_ranger.tif'
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
                  useiso3 = 'IDN', 
                  admin_unit_level = 'ADMIN2',
                  pr_country = 'country',
                  api_year = 2012)

data_ml_cov <- load_data(PR_path, 
                         API_path, 
                         pop_path, 
                         ml_local_raster_paths, 
                         shapefile_path, 
                         shapefile_pattern = '.shp$', 
                         useiso3 = 'IDN', 
                         admin_unit_level = 'ADMIN2',
                         pr_country = 'country',
                         api_year = 2012)

data_all_cov <- load_data(PR_path, 
                          API_path, 
                          pop_path, 
                          c(cov_raster_paths, ml_local_raster_paths), 
                          shapefile_path, 
                          shapefile_pattern = '.shp$', 
                          useiso3 = 'IDN', 
                          admin_unit_level = 'ADMIN2',
                          pr_country = 'country',
                          api_year = 2012)


# pre analysis

data_idn_cov <- process_data(
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
  useiso3 = 'IDN',
  transform = c(4:7))

save(data_idn_cov, file = 'model_outputs/idn_cov_data.RData')

data_idn_ml <- process_data(
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
  useiso3 = 'IDN',
  transform = NULL)
save(data_idn_ml, file = 'model_outputs/idn_ml_data.RData')

data_idn_all <- process_data(
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
  useiso3 = 'IDN',
  transform = c(4:7))

save(data_idn_all, file = 'model_outputs/idn_all_data.RData')


autoplot(data_idn_cov, pr_limits = c(0, 0.3))
autoplot(data_idn_cov, pr_limits = c(0, 0.3), trans = 'log1p')

ggsave('figs/idn_input_data.png')



mesh_idn <- build_mesh(data_idn_cov, mesh.args = list(max.edge = c(0.7, 5), cut = 0.7))
autoplot(mesh_idn)
save(mesh_idn, file = 'model_outputs/idn_mesh.RData')



# Define cross validation strategies
data_cv1_idn <- cv_random_folds(data_idn_cov, k = 6)
data_cv1_idn_ml <- cv_random_folds(data_idn_ml, k = 6, 
                                   polygon_folds = attr(data_cv1_idn, 'polygon_folds'),
                                   pr_folds = attr(data_cv1_idn, 'pr_folds'))
data_cv1_idn_all <- cv_random_folds(data_idn_all, k = 6, 
                                    polygon_folds = attr(data_cv1_idn, 'polygon_folds'),
                                    pr_folds = attr(data_cv1_idn, 'pr_folds'))


autoplot(data_cv1_idn, jitter = 0)
autoplot(data_cv1_idn_ml, jitter = 0)

ggsave('figs/idn_cv_random.png')
save(data_cv1_idn, file = 'model_outputs/idn_cv_1.RData')
save(data_cv1_idn_ml, file = 'model_outputs/idn_cv_1_ml.RData')
save(data_cv1_idn_all, file = 'model_outputs/idn_cv_1_all.RData')


# Spatial
data_cv2_idn <- cv_spatial_folds(data_idn_cov, k = 6)
data_cv2_idn_ml <- cv_spatial_folds(data_idn_ml, k = 6, 
                                    polygon_folds = attr(data_cv2_idn, 'polygon_folds'),
                                    pr_folds = attr(data_cv2_idn, 'pr_folds'))
data_cv2_idn_all <- cv_spatial_folds(data_idn_all, k = 6, 
                                     polygon_folds = attr(data_cv2_idn, 'polygon_folds'),
                                     pr_folds = attr(data_cv2_idn, 'pr_folds'))

autoplot(data_cv2_idn, jitter = 0)
autoplot(data_cv2_idn_ml, jitter = 0)

ggsave('figs/idn_cv_spatial2.png')
save(data_cv2_idn, file = 'model_outputs/idn_cv_2.RData')
save(data_cv2_idn_ml, file = 'model_outputs/idn_cv_2_ml.RData')
save(data_cv2_idn_all, file = 'model_outputs/idn_cv_2_all.RData')

#autoplot(data_cv1_idn[[1]]$train, pr_limits = c(0, 0.3))






