



run_machine_learn <- function(pr_path, 
                              covs, 
                              pr_region, 
                              pr_min_year = 1990,
                              extent,
                              model_list, 
                              search_vec,
                              tuneLength_vec,
                              pr_pos_column = 'positive',
                              pr_n_column = 'examined',
                              pr_latlon = c('latitude', 'longitude'),
                              pr_country = 'country',
                              pr_age_low = 'lower_age',
                              pr_age_high = 'upper_age',
                              standardisePR = c(2, 10),
                              roundPR = FALSE,
                              standardisePars = 'Pf_Smith2007',
                              ){
  
  
  
  pr <- readr::read_csv(PR_path, guess_max  = 1e5)
  
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
  
  models <- list()

  y <- pr_clean$prevalence
  partition <- createMultiFolds(y, k = 5, times = 1)
  


  models[[1]] <- train(pr_extracted, y, 
                       method = 'enet',
                       trControl = trainControl(index = partition, 
                                                returnData = TRUE,
                                                savePredictions = TRUE, 
                                                search = search_vec[1]),
                       tuneLength = tuneLength_vec[1])
                       



  models[[2]] <- train(pr_extracted, y, 
                       method = 'xgbTree',
                       trControl = trainControl(index = partition, 
                                                returnData = TRUE,
                                                savePredictions = TRUE, 
                                                search = search_vec[2]),
                       tuneLength = tuneLength_vec[2])
                       





  models[[3]] <- train(pr_extracted, y, 
                       method = 'ranger',
                       trControl = trainControl(index = partition, 
                                                returnData = TRUE,
                                                savePredictions = TRUE, 
                                                search = search_vec[3]),
                       tuneLength = tuneLength_vec[3])
                       


  models[[4]] <- train(pr_extracted, y, 
                       method = 'ppr',
                       trControl = trainControl(index = partition, 
                                                returnData = TRUE,
                                                savePredictions = TRUE, 
                                                search = search_vec[4]),
                       tuneLength = tuneLength_vec[4])
                       

  
  if(!is.null(figpath)){
    png(paste0(figpath, 'enetopt.png'))
    print(plot(models[[1]]))
    dev.off()
  } 
  

    print(plotObsVsPred(extractPrediction(list(models[[1]]))))
    print(plotObsVsPred(extractPrediction(list(models[[2]]))))
    print(plotObsVsPred(extractPrediction(list(models[[3]]))))
    print(plotObsVsPred(extractPrediction(list(models[[4]]))))


  covs_crop <- crop(covs, extent)


  return(machine_learn_covariates)
  
}
















  
