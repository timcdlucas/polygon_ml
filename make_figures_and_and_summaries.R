
# figure 1.cross validation. %% Do fig 1 and 2, random and spatial cv. IDN on top, MDG and SEN below in each.
# figure 3 and 4. data and predicted incidence maps. Indonesia and Senegal only. Fig 3 ind, fig 4 sen Data, Rand, Spatial for best model? Joint model?
# figure 5, 6. Spat and random cv. PR vs Poly columns, countries as rows, model as colour?


if(Sys.info()["user"] != 'anita'){
  setwd('~/timz/timothy/polygon_ml_wsc')
} else {
  setwd('~/Z/timothy/polygon_ml_wsc')
}

# Libs

## Spatial packages
library(raster)
library(maptools)
library(rgeos)

## dataframe packages
library(dplyr)
library(readr)
library(magrittr)
library(tidyr)

library(malariaAtlas)

library(ggplot2)
library(cowplot)
library(scales) 
theme_set(theme_minimal())
#theme_update(text = element_text(size = 10))


# Metrics 



source('plotting_functions.R')

# Paths


## MDG


### Cross validation object

data_cv1_mdg_path <- 'model_outputs/mdg_cv_1.RData'
data_cv2_mdg_path <- 'model_outputs/mdg_cv_2.RData'


### CV 1 output

cv1_covs_mdg_path <- 'model_outputs/mdg_covs_cv_1.RData'
cv1_ml_mdg_path <- 'model_outputs/mdg_ml_cv_1.RData'
cv1_all_mdg_path <- 'model_outputs/mdg_all_cv_1.RData'


### CV 2 output

cv2_covs_mdg_path <- 'model_outputs/mdg_covs_cv_2.RData'
cv2_ml_mdg_path <- 'model_outputs/mdg_ml_cv_2.RData'
cv2_all_mdg_path <- 'model_outputs/mdg_all_cv_2.RData'





### Cross validation object

data_cv1_col_path <- 'model_outputs/col_cv_1.RData'
data_cv2_col_path <- 'model_outputs/col_cv_2.RData'


### CV 1 output

cv1_covs_col_path <- 'model_outputs/col_covs_cv_1.RData'
cv1_ml_col_path <- 'model_outputs/col_ml_cv_1.RData'
cv1_all_col_path <- 'model_outputs/col_all_cv_1.RData'


### CV 2 output

cv2_covs_col_path <- 'model_outputs/col_covs_cv_2.RData'
cv2_ml_col_path <- 'model_outputs/col_ml_cv_2.RData'
cv2_all_col_path <- 'model_outputs/col_all_cv_2.RData'




# Read all data



### Cross validation object

load(data_cv1_mdg_path)
load(data_cv2_mdg_path)


### CV 1 output

cv1_covs_mdg <- get(load(cv1_covs_mdg_path))
cv1_ml_mdg <- get(load(cv1_ml_mdg_path))
cv1_all_mdg <- get(load(cv1_all_mdg_path))


### CV 2 output

cv2_covs_mdg <- get(load(cv2_covs_mdg_path))
cv2_ml_mdg <- get(load(cv2_ml_mdg_path))
cv2_all_mdg <- get(load(cv2_all_mdg_path))





### Cross validation object

load(data_cv1_col_path)
load(data_cv2_col_path)


### CV 1 output

cv1_covs_col <- get(load(cv1_covs_col_path))
cv1_ml_col <- get(load(cv1_ml_col_path))
cv1_all_col <- get(load(cv1_all_col_path))


### CV 2 output

cv2_covs_col <- get(load(cv2_covs_col_path))
cv2_ml_col <- get(load(cv2_ml_col_path))
cv2_all_col <- get(load(cv2_all_col_path))









# Join all table data.

cv1_df_mdg <- rbind(cv1_covs_mdg$summary$combined_aggregated %>% mutate(model = 'Covariates only'),
                    cv1_ml_mdg$summary$combined_aggregated %>% mutate(model = 'ML models only'),
                    cv1_all_mdg$summary$combined_aggregated %>% mutate(model = 'Both')) %>% 
                 mutate(country = 'Madagascar')

cv2_df_mdg <- rbind(cv2_covs_mdg$summary$combined_aggregated %>% mutate(model = 'Covariates only'),
                    cv2_ml_mdg$summary$combined_aggregated %>% mutate(model = 'ML models only'),
                    cv2_all_mdg$summary$combined_aggregated %>% mutate(model = 'Both')) %>% 
                 mutate(country = 'Madagascar')

cv1_df_col <- rbind(cv1_covs_col$summary$combined_aggregated %>% mutate(model = 'Covariates only'),
                    cv1_ml_col$summary$combined_aggregated %>% mutate(model = 'ML models only'),
                    cv1_all_col$summary$combined_aggregated %>% mutate(model = 'Both')) %>% 
                 mutate(country = 'Colombia')

cv2_df_col <- rbind(cv2_covs_col$summary$combined_aggregated %>% mutate(model = 'Covariates only'),
                    cv2_ml_col$summary$combined_aggregated %>% mutate(model = 'ML models only'),
                    cv2_all_col$summary$combined_aggregated %>% mutate(model = 'Both')) %>% 
                 mutate(country = 'Colombia')


cv1_df <- rbind(cv1_df_mdg, cv1_df_col) %>% 
  mutate(model = factor(model, levels = c('Covariates only', 'ML models only', 'Both')))
cv2_df <- rbind(cv2_df_mdg, cv2_df_col) %>% 
  mutate(model = factor(model, levels = c('Covariates only', 'ML models only', 'Both')))






# CV plots


mysqrt_trans <- function() { 
   trans_new("mysqrt", 
             transform = base::sqrt, 
             inverse = function(x) ifelse(x<0, 0, x^2), 
             domain = c(0, Inf)) 
} 


  p1 <- 
    ggplot(cv1_df, aes(response, pred_api)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = 'grey', size = 4) +
      geom_point(size = 6, alpha = 0.6) +
      facet_grid(country ~ model) +
      scale_x_continuous(trans = "mysqrt", breaks = c(0, 20, 100, 300, 500)) + 
      scale_y_continuous(trans = "mysqrt", breaks = c(0, 5, 20, 50, 100, 200, 300)) +
      xlab('Observed API') + 
      ylab('Predicted API') + 
      theme_bw() +
      theme(text = element_text(size = 50),
            panel.spacing = unit(3, "lines"))

  png('figs/cv1_scatter.png', width = 1600, height = 1100)
  print(p1)
  dev.off()




p2 <- 
  ggplot(cv2_df, aes(response, pred_api, colour = factor(fold))) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = 'grey', size = 4) +
      geom_point(size = 6, alpha = 0.6) +
      facet_grid(country ~ model, scale = 'free_y') +
      scale_x_continuous(trans = "mysqrt", breaks = c(0, 20, 100, 300, 500)) + 
      scale_y_continuous(trans = "mysqrt", breaks = c(0, 5, 20, 50, 100, 200, 300)) +
      xlab('Observed API') + 
      ylab('Predicted API') + 
      theme_bw() +
      guides(colour = FALSE) +
      theme(text = element_text(size = 50),
            panel.spacing = unit(3, "lines"))

png('figs/cv2_scatter.png', width = 1600, height = 1100)
print(p2)
dev.off()





p3 <- obspred_map(data_cv1_col, cv1_ml_col, trans = 'log10', column = FALSE, mask = TRUE)

p4 <- list()

p4[[1]] <- 
  p3[[1]] +
    guides(fill = FALSE) + 
    xlab('Longitude') + 
    ylab('Latitude') +
    ggtitle('Observed data') +
    theme(text = element_text(size = 30)) +
    xlim(-80, NA) + 
    ylim(-4.5, 13)



p4[[2]] <- 
  p3[[2]] +
    xlab('Longitude') + 
    ylab('') +
    theme(text = element_text(size = 30)) +
    ggtitle('Out-of-sample predictions') +
    scale_fill_viridis_c(trans = 'log10', 
                         breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3),
                         labels = trans_format("log10", math_format(10^.x))) +
    guides(fill = guide_colourbar(barheight = 30, title="API")) 


png('figs/col_obs_pred_map_ml.png', width = 1600, height = 1200)
print(plot_grid(plotlist = p4, rel_widths = c(1, 1.07)))
dev.off()






