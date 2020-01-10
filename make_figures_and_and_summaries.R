
# figure 1.cross validation. %% Do fig 1 and 2, random and spatial cv. IDN on top, MDG and SEN below in each.
# figure 3 and 4. data and predicted incidence maps. Indonesia and Senegal only. Fig 3 ind, fig 4 sen Data, Rand, Spatial for best model? Joint model?
# figure 5, 6. Spat and random cv. PR vs Poly columns, countries as rows, model as colour?


if(Sys.info()["user"] != 'anita'){
  setwd('~/timz/timothy/polygon_ml')
} else {
  setwd('~/Z/timothy/polygon_ml')
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


mysqrt_trans <- function() { 
  trans_new("mysqrt", 
            transform = base::sqrt, 
            inverse = function(x) ifelse(x<0, 0, x^2), 
            domain = c(0, Inf)) 
} 

# Metrics 



source('plotting_functions.R')

# Paths

# --------------------------------------------------------------------- #
## MDG


### Cross validation object

data_cv1_mdg_path <- 'model_outputs/mdg_cv_1.RData'
data_cv2_mdg_path <- 'model_outputs/mdg_cv_2.RData'


### CV 1 output

cv1_covs_mdg_path <- 'model_outputs/mdg_covs_cv_1.RData'
cv1_ml_mdg_path <- 'model_outputs/mdg_ml_cv_1.RData'
cv1_all_mdg_path <- 'model_outputs/mdg_all_cv_1.RData'
cv1_mle_mdg_path <- 'model_outputs/mdg_mle_cv_1.RData'
cv1_mlg_mdg_path <- 'model_outputs/mdg_mlg_cv_1.RData'


### CV 2 output

cv2_covs_mdg_path <- 'model_outputs/mdg_covs_cv_2.RData'
cv2_ml_mdg_path <- 'model_outputs/mdg_ml_cv_2.RData'
cv2_all_mdg_path <- 'model_outputs/mdg_all_cv_2.RData'
cv2_mle_mdg_path <- 'model_outputs/mdg_mle_cv_2.RData'
cv2_mlg_mdg_path <- 'model_outputs/mdg_mlg_cv_2.RData'

# --------------------------------------------------------------------- #


# --------------------------------------------------------------------- #
# col

### Cross validation object

data_cv1_col_path <- 'model_outputs/col_cv_1.RData'
data_cv2_col_path <- 'model_outputs/col_cv_2.RData'


### CV 1 output

cv1_covs_col_path <- 'model_outputs/col_covs_cv_1.RData'
cv1_ml_col_path <- 'model_outputs/col_ml_cv_1.RData'
cv1_all_col_path <- 'model_outputs/col_all_cv_1.RData'
cv1_mle_col_path <- 'model_outputs/col_mle_cv_1.RData'
cv1_mlg_col_path <- 'model_outputs/col_mlg_cv_1.RData'

### CV 2 output

cv2_covs_col_path <- 'model_outputs/col_covs_cv_2.RData'
cv2_ml_col_path <- 'model_outputs/col_ml_cv_2.RData'
cv2_all_col_path <- 'model_outputs/col_all_cv_2.RData'
cv2_mle_col_path <- 'model_outputs/col_mle_cv_2.RData'
cv2_mlg_col_path <- 'model_outputs/col_mlg_cv_2.RData'
# --------------------------------------------------------------------- #


# --------------------------------------------------------------------- #
# sen

### Cross validation object

data_cv1_sen_path <- 'model_outputs/sen_cv_1.RData'
data_cv2_sen_path <- 'model_outputs/sen_cv_2.RData'


### CV 1 output

cv1_covs_sen_path <- 'model_outputs/sen_covs_cv_1.RData'
cv1_ml_sen_path <- 'model_outputs/sen_ml_cv_1.RData'
cv1_all_sen_path <- 'model_outputs/sen_all_cv_1.RData'
cv1_mle_sen_path <- 'model_outputs/sen_mle_cv_1.RData'
cv1_mlg_sen_path <- 'model_outputs/sen_mlg_cv_1.RData'

### CV 2 output

cv2_covs_sen_path <- 'model_outputs/sen_covs_cv_2.RData'
cv2_ml_sen_path <- 'model_outputs/sen_ml_cv_2.RData'
cv2_all_sen_path <- 'model_outputs/sen_all_cv_2.RData'
cv2_mle_sen_path <- 'model_outputs/sen_mle_cv_2.RData'
cv2_mlg_sen_path <- 'model_outputs/sen_mlg_cv_2.RData'
# --------------------------------------------------------------------- #



# --------------------------------------------------------------------- #
# idn

### Cross validation object

data_cv1_idn_path <- 'model_outputs/idn_cv_1.RData'
data_cv2_idn_path <- 'model_outputs/idn_cv_2.RData'


### CV 1 output

cv1_covs_idn_path <- 'model_outputs/idn_covs_cv_1.RData'
cv1_ml_idn_path <- 'model_outputs/idn_ml_cv_1.RData'
cv1_all_idn_path <- 'model_outputs/idn_all_cv_1.RData'
cv1_mle_idn_path <- 'model_outputs/idn_mle_cv_1.RData'
cv1_mlg_idn_path <- 'model_outputs/idn_mlg_cv_1.RData'

### CV 2 output

cv2_covs_idn_path <- 'model_outputs/idn_covs_cv_2.RData'
cv2_ml_idn_path <- 'model_outputs/idn_ml_cv_2.RData'
cv2_all_idn_path <- 'model_outputs/idn_all_cv_2.RData'
cv2_mle_idn_path <- 'model_outputs/idn_mle_cv_2.RData'
cv2_mlg_idn_path <- 'model_outputs/idn_mlg_cv_2.RData'
# --------------------------------------------------------------------- #



# Read all data



### Cross validation object

load(data_cv1_mdg_path)
load(data_cv2_mdg_path)


### CV 1 output

cv1_covs_mdg <- get(load(cv1_covs_mdg_path))
cv1_ml_mdg <- get(load(cv1_ml_mdg_path))
cv1_all_mdg <- get(load(cv1_all_mdg_path))
cv1_mle_mdg <- get(load(cv1_mle_mdg_path))
cv1_mlg_mdg <- get(load(cv1_mlg_mdg_path))

### CV 2 output

cv2_covs_mdg <- get(load(cv2_covs_mdg_path))
cv2_ml_mdg <- get(load(cv2_ml_mdg_path))
cv2_all_mdg <- get(load(cv2_all_mdg_path))
cv2_mle_mdg <- get(load(cv2_mle_mdg_path))
cv2_mlg_mdg <- get(load(cv2_mlg_mdg_path))



### Cross validation object

load(data_cv1_col_path)
load(data_cv2_col_path)


### CV 1 output

cv1_covs_col <- get(load(cv1_covs_col_path))
cv1_ml_col <- get(load(cv1_ml_col_path))
cv1_all_col <- get(load(cv1_all_col_path))
cv1_mle_col <- get(load(cv1_mle_col_path))
cv1_mlg_col <- get(load(cv1_mlg_col_path))

### CV 2 output

cv2_covs_col <- get(load(cv2_covs_col_path))
cv2_ml_col <- get(load(cv2_ml_col_path))
cv2_all_col <- get(load(cv2_all_col_path))
cv2_mle_col <- get(load(cv2_mle_col_path))
cv2_mlg_col <- get(load(cv2_mlg_col_path))




### Cross validation object

load(data_cv1_sen_path)
load(data_cv2_sen_path)


### CV 1 output

cv1_covs_sen <- get(load(cv1_covs_sen_path))
cv1_ml_sen <- get(load(cv1_ml_sen_path))
cv1_all_sen <- get(load(cv1_all_sen_path))
cv1_mle_sen <- get(load(cv1_mle_sen_path))
cv1_mlg_sen <- get(load(cv1_mlg_sen_path))

### CV 2 output

cv2_covs_sen <- get(load(cv2_covs_sen_path))
cv2_ml_sen <- get(load(cv2_ml_sen_path))
cv2_all_sen <- get(load(cv2_all_sen_path))
cv2_mle_sen <- get(load(cv2_mle_sen_path))
cv2_mlg_sen <- get(load(cv2_mlg_sen_path))




### Cross validation object

load(data_cv1_idn_path)
load(data_cv2_idn_path)


### CV 1 output

cv1_covs_idn <- get(load(cv1_covs_idn_path))
cv1_ml_idn <- get(load(cv1_ml_idn_path))
cv1_all_idn <- get(load(cv1_all_idn_path))
cv1_mle_idn <- get(load(cv1_mle_idn_path))
cv1_mlg_idn <- get(load(cv1_mlg_idn_path))

### CV 2 output

cv2_covs_idn <- get(load(cv2_covs_idn_path))
cv2_ml_idn <- get(load(cv2_ml_idn_path))
cv2_all_idn <- get(load(cv2_all_idn_path))
cv2_mle_idn <- get(load(cv2_mle_idn_path))
cv2_mlg_idn <- get(load(cv2_mlg_idn_path))







# Join all table data.

cv1_df_mdg <- rbind(cv1_covs_mdg$summary$combined_aggregated %>% mutate(model = 'Envir'),
                    cv1_ml_mdg$summary$combined_aggregated %>% mutate(model = 'MLl'),
                    cv1_mlg_mdg$summary$combined_aggregated %>% mutate(model = 'MLg'),
                    cv1_all_mdg$summary$combined_aggregated %>% mutate(model = 'MLl + MLg'),
                    cv1_mle_mdg$summary$combined_aggregated %>% mutate(model = 'Envir + MLl')) %>% 
                 mutate(country = 'MDG')

cv2_df_mdg <- rbind(cv2_covs_mdg$summary$combined_aggregated %>% mutate(model = 'Envir'),
                    cv2_ml_mdg$summary$combined_aggregated %>% mutate(model = 'MLl'),
                    cv2_mlg_mdg$summary$combined_aggregated %>% mutate(model = 'MLg'),
                    cv2_all_mdg$summary$combined_aggregated %>% mutate(model = 'MLl + MLg'),
                    cv2_mle_mdg$summary$combined_aggregated %>% mutate(model = 'Envir + MLl')) %>% 
                 mutate(country = 'MDG')


cv1_df_col <- rbind(cv1_covs_col$summary$combined_aggregated %>% mutate(model = 'Envir'),
                    cv1_ml_col$summary$combined_aggregated %>% mutate(model = 'MLl'),
                    cv1_mlg_col$summary$combined_aggregated %>% mutate(model = 'MLg'),
                    cv1_all_col$summary$combined_aggregated %>% mutate(model = 'MLl + MLg'),
                    cv1_mle_col$summary$combined_aggregated %>% mutate(model = 'Envir + MLl')) %>% 
                 mutate(country = 'COL')

cv2_df_col <- rbind(cv2_covs_col$summary$combined_aggregated %>% mutate(model = 'Envir'),
                    cv2_ml_col$summary$combined_aggregated %>% mutate(model = 'MLl'),
                    cv2_mlg_col$summary$combined_aggregated %>% mutate(model = 'MLg'),
                    cv2_all_col$summary$combined_aggregated %>% mutate(model = 'MLl + MLg'),
                    cv2_mle_col$summary$combined_aggregated %>% mutate(model = 'Envir + MLl')) %>% 
                 mutate(country = 'COL')


cv1_df_sen <- rbind(cv1_covs_sen$summary$combined_aggregated %>% mutate(model = 'Envir'),
                    cv1_ml_sen$summary$combined_aggregated %>% mutate(model = 'MLl'),
                    cv1_mlg_sen$summary$combined_aggregated %>% mutate(model = 'MLg'),
                    cv1_all_sen$summary$combined_aggregated %>% mutate(model = 'MLl + MLg'),
                    cv1_mle_sen$summary$combined_aggregated %>% mutate(model = 'Envir + MLl')) %>% 
                 mutate(country = 'SEN')

cv2_df_sen <- rbind(cv2_covs_sen$summary$combined_aggregated %>% mutate(model = 'Envir'),
                    cv2_ml_sen$summary$combined_aggregated %>% mutate(model = 'MLl'),
                    cv2_mlg_sen$summary$combined_aggregated %>% mutate(model = 'MLg'),
                    cv2_all_sen$summary$combined_aggregated %>% mutate(model = 'MLl + MLg'),
                    cv2_mle_sen$summary$combined_aggregated %>% mutate(model = 'Envir + MLl')) %>% 
                 mutate(country = 'SEN')


cv1_df_idn <- rbind(cv1_covs_idn$summary$combined_aggregated %>% mutate(model = 'Envir'),
                    cv1_ml_idn$summary$combined_aggregated %>% mutate(model = 'MLl'),
                    cv1_mlg_idn$summary$combined_aggregated %>% mutate(model = 'MLg'),
                    cv1_all_idn$summary$combined_aggregated %>% mutate(model = 'MLl + MLg'),
                    cv1_mle_idn$summary$combined_aggregated %>% mutate(model = 'Envir + MLl')) %>% 
                 mutate(country = 'IDN')

cv2_df_idn <- rbind(cv2_covs_idn$summary$combined_aggregated %>% mutate(model = 'Envir'),
                    cv2_ml_idn$summary$combined_aggregated %>% mutate(model = 'MLl'),
                    cv2_mlg_idn$summary$combined_aggregated %>% mutate(model = 'MLg'),
                    cv2_all_idn$summary$combined_aggregated %>% mutate(model = 'MLl + MLg'),
                    cv2_mle_idn$summary$combined_aggregated %>% mutate(model = 'Envir + MLl')) %>% 
                 mutate(country = 'IDN')




cv1_df <- rbind(cv1_df_mdg, cv1_df_col, cv1_df_sen, cv1_df_idn) %>% 
  mutate(model = factor(model, levels = c('Envir', 'MLl', 'Envir + MLl', 'MLg', 'MLl + MLg')))
cv2_df <- rbind(cv2_df_mdg, cv2_df_col, cv2_df_sen, cv2_df_idn) %>% 
  mutate(model = factor(model, levels = c('Envir', 'MLl', 'Envir + MLl', 'MLg', 'MLl + MLg')))


write.csv(cv1_df, 'figs/random_all_estiamtes.csv')
write.csv(cv2_df, 'figs/spatial_all_estiamtes.csv')


# CV plots



  p1 <- 
    cv1_df %>% 
      mutate(model_country = paste(country, model)) %>% 
    ggplot(aes(response, pred_api)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = 'grey', size = 4) +
      geom_point(size = 6, alpha = 0.6) +
      facet_wrap(model_country ~ ., scales = 'free') +
      scale_x_continuous(trans = "mysqrt", breaks = c(0, 20, 100, 300, 500)) + 
      scale_y_continuous(trans = "mysqrt", breaks = c(0, 5, 20, 50, 100, 200, 300)) +
      xlab('Observed API') + 
      ylab('Predicted API') + 
      theme_bw() +
      theme(text = element_text(size = 40),
            panel.spacing = unit(1, "lines"))

  png('figs/cv1_scatter.png', width = 1900, height = 1500)
  print(p1)
  dev.off()




p2 <- 
  cv2_df %>% 
    mutate(model_country = paste(country, model)) %>% 
  ggplot(aes(response, pred_api)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = 'grey', size = 4) +
      geom_point(size = 6, alpha = 0.6) +
      facet_wrap(model_country ~ ., scales = 'free') +
      scale_x_continuous(trans = "mysqrt", breaks = c(0, 20, 100, 300, 500)) + 
      scale_y_continuous(trans = "mysqrt", breaks = c(0, 5, 20, 50, 100, 200, 300)) +
      xlab('Observed API') + 
      ylab('Predicted API') + 
      theme_bw() +
      guides(colour = FALSE) +
      theme(text = element_text(size = 40),
            panel.spacing = unit(1, "lines"))

png('figs/cv2_scatter.png', width = 1900, height = 1500)
print(p2)
dev.off()



## CV plots with less per plot.

cv1_df_l <- cv1_df %>% filter(model %in% c('Envir', 'MLl', 'Envir + MLl'))
cv2_df_l <- cv2_df %>% filter(model %in% c('Envir', 'MLl', 'Envir + MLl'))

res1 <-
  c(cv1_covs_mdg$summary$polygon_metrics$pearson,
    cv1_ml_mdg$summary$polygon_metrics$pearson,
    cv1_mle_mdg$summary$polygon_metrics$pearson,
    cv1_covs_col$summary$polygon_metrics$pearson,
    cv1_ml_col$summary$polygon_metrics$pearson,
    cv1_mle_col$summary$polygon_metrics$pearson,
    cv1_covs_sen$summary$polygon_metrics$pearson,
    cv1_ml_sen$summary$polygon_metrics$pearson,
    cv1_mle_sen$summary$polygon_metrics$pearson,
    cv1_covs_idn$summary$polygon_metrics$pearson,
    cv1_ml_idn$summary$polygon_metrics$pearson,
    cv1_mle_idn$summary$polygon_metrics$pearson)

cv1_l_labels <- 
  cv1_df_l %>% 
    distinct(country, model) %>% 
    #mutate(model_country = paste0(country, ': ',  model)) %>% 
    mutate(corr = as.character(sprintf("%.2f", round(res1, 2))))

spread(cv1_l_labels, model, corr)

cv1_df_l <- left_join(cv1_df_l, cv1_l_labels)


p1 <- 
  cv1_df_l %>% 
  mutate(model_country = paste0(country, ', ',  model, ': ', corr)) %>% 
  ggplot(aes(response, pred_api)) +
  geom_point(size = 6, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 5, colour = 'grey', size = 3) +
  facet_wrap(model_country ~ ., scales = 'free', ncol = 3) +
  scale_x_continuous(trans = "mysqrt", breaks = c(0, 20, 100, 200, 400)) + 
  scale_y_continuous(trans = "mysqrt", breaks = c(0, 5, 20, 50, 100, 200, 400)) +
  xlab('Observed API') + 
  ylab('Predicted API') + 
  theme_bw() +
  theme(text = element_text(size = 40),
        panel.spacing = unit(1, "lines"))

png('figs/cv1_l_scatter.png', width = 1300, height = 1500)
print(p1)
dev.off()



res2 <-
  c(cv2_covs_mdg$summary$polygon_metrics$pearson,
    cv2_ml_mdg$summary$polygon_metrics$pearson,
    cv2_mle_mdg$summary$polygon_metrics$pearson,
    cv2_covs_col$summary$polygon_metrics$pearson,
    cv2_ml_col$summary$polygon_metrics$pearson,
    cv2_mle_col$summary$polygon_metrics$pearson,
    cv2_covs_sen$summary$polygon_metrics$pearson,
    cv2_ml_sen$summary$polygon_metrics$pearson,
    cv2_mle_sen$summary$polygon_metrics$pearson,
    cv2_covs_idn$summary$polygon_metrics$pearson,
    cv2_ml_idn$summary$polygon_metrics$pearson,
    cv2_mle_idn$summary$polygon_metrics$pearson)

cv2_l_labels <- 
  cv2_df_l %>% 
  distinct(country, model) %>% 
  #mutate(model_country = paste0(country, ': ',  model)) %>% 
  mutate(corr = as.character(sprintf("%.2f", round(res2, 2))))

spread(cv2_l_labels, model, corr)


cv2_df_l <- left_join(cv2_df_l, cv2_l_labels)


p2 <- 
  cv2_df_l %>% 
  mutate(model_country = paste0(country, ', ',  model, ': ', corr)) %>% 
  ggplot(aes(response, pred_api)) +
  geom_point(size = 6, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 5, colour = 'grey', size = 3) +
  facet_wrap(model_country ~ ., scales = 'free', ncol = 3) +
  scale_x_continuous(trans = "mysqrt", breaks = c(0, 20, 100, 200, 400)) + 
  scale_y_continuous(trans = "mysqrt", breaks = c(0, 5, 20, 50, 100, 200, 400)) +
  xlab('Observed API') + 
  ylab('Predicted API') + 
  theme_bw() +
  guides(colour = FALSE) +
  theme(text = element_text(size = 40),
        panel.spacing = unit(1, "lines"))

png('figs/cv2_l_scatter.png', width = 1300, height = 1500)
print(p2)
dev.off()






cv1_df_g <- cv1_df %>% filter(model %in% c('Envir', 'MLg', 'MLl + MLg'))
cv2_df_g <- cv2_df %>% filter(model %in% c('Envir', 'MLg', 'MLl + MLg'))



p1 <- 
  cv1_df_g %>% 
  mutate(model_country = paste(country, model)) %>% 
  ggplot(aes(response, pred_api)) +
  geom_point(size = 6, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 5, colour = 'grey', size = 3) +
  facet_wrap(model_country ~ ., scales = 'free', ncol = 3) +
  scale_x_continuous(trans = "mysqrt", breaks = c(0, 20, 100, 200, 400)) + 
  scale_y_continuous(trans = "mysqrt", breaks = c(0, 5, 20, 50, 100, 200, 400)) +
  xlab('Observed API') + 
  ylab('Predicted API') + 
  theme_bw() +
  theme(text = element_text(size = 40),
        panel.spacing = unit(1, "lines"))

png('figs/cv1_g_scatter.png', width = 1300, height = 1500)
print(p1)
dev.off()




p2 <- 
  cv2_df_g %>% 
  mutate(model_country = paste(country, model)) %>% 
  ggplot(aes(response, pred_api)) +
  geom_point(size = 6, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 5, colour = 'grey', size = 3) +
  facet_wrap(model_country ~ ., scales = 'free', ncol = 3) +
  scale_x_continuous(trans = "mysqrt", breaks = c(0, 20, 100, 200, 400)) + 
  scale_y_continuous(trans = "mysqrt", breaks = c(0, 5, 20, 50, 100, 200, 400)) +
  xlab('Observed API') + 
  ylab('Predicted API') + 
  theme_bw() +
  guides(colour = FALSE) +
  theme(text = element_text(size = 40),
        panel.spacing = unit(1, "lines"))

png('figs/cv2_g_scatter.png', width = 1300, height = 1500)
print(p2)
dev.off()



# ----------     #
# global results # 
# -------------- #

res3 <-
  c(cv1_covs_mdg$summary$polygon_metrics$pearson,
    cv1_mlg_mdg$summary$polygon_metrics$pearson,
    cv1_all_mdg$summary$polygon_metrics$pearson,
    cv1_covs_col$summary$polygon_metrics$pearson,
    cv1_mlg_col$summary$polygon_metrics$pearson,
    cv1_all_col$summary$polygon_metrics$pearson,
    cv1_covs_sen$summary$polygon_metrics$pearson,
    cv1_mlg_sen$summary$polygon_metrics$pearson,
    cv1_all_sen$summary$polygon_metrics$pearson,
    cv1_covs_idn$summary$polygon_metrics$pearson,
    cv1_mlg_idn$summary$polygon_metrics$pearson,
    cv1_all_idn$summary$polygon_metrics$pearson)


cv1_g_labels <- 
  cv1_df_g %>% 
  distinct(country, model) %>% 
  #mutate(model_country = paste0(country, ': ',  model)) %>% 
  mutate(corr = as.character(sprintf("%.2f", round(res3, 2))))

spread(cv1_g_labels, model, corr)

res4 <-
  c(cv2_covs_mdg$summary$polygon_metrics$pearson,
    cv2_mlg_mdg$summary$polygon_metrics$pearson,
    cv2_all_mdg$summary$polygon_metrics$pearson,
    cv2_covs_col$summary$polygon_metrics$pearson,
    cv2_mlg_col$summary$polygon_metrics$pearson,
    cv2_all_col$summary$polygon_metrics$pearson,
    cv2_covs_sen$summary$polygon_metrics$pearson,
    cv2_mlg_sen$summary$polygon_metrics$pearson,
    cv2_all_sen$summary$polygon_metrics$pearson,
    cv2_covs_idn$summary$polygon_metrics$pearson,
    cv2_mlg_idn$summary$polygon_metrics$pearson,
    cv2_all_idn$summary$polygon_metrics$pearson)



cv2_g_labels <- 
  cv2_df_g %>% 
  distinct(country, model) %>% 
  #mutate(model_country = paste0(country, ': ',  model)) %>% 
  mutate(corr = as.character(sprintf("%.2f", round(res4, 2))))

spread(cv2_g_labels, model, corr)


##---------------------------------------##
## Maps                                  ##  
##---------------------------------------##

p3 <- obspred_map(data_cv2_col, cv2_ml_col, trans = 'log10', column = FALSE, mask = TRUE)

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



p3 <- obspred_map(data_cv1_mdg, cv1_ml_mdg, trans = 'log10', column = FALSE, mask = TRUE)

p4 <- list()

p4[[1]] <- 
  p3[[1]] +
    guides(fill = FALSE) + 
    xlab('Longitude') + 
    ylab('Latitude') +
    ggtitle('Observed data') +
    theme(text = element_text(size = 30)) 



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


png('figs/mdg_obs_pred_map_ml.png', width = 1600, height = 1200)
print(plot_grid(plotlist = p4, rel_widths = c(1, 1.07)))
dev.off()



p3 <- obspred_map(data_cv1_idn, cv1_ml_idn, trans = 'log10', column = FALSE, mask = TRUE)

p4 <- list()

p4[[1]] <- 
  p3[[1]] +
    guides(fill = FALSE) + 
    xlab('Longitude') + 
    ylab('Latitude') +
    ggtitle('Observed data') +
    theme(text = element_text(size = 30)) 


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


png('figs/idn_obs_pred_map_ml.png', width = 1600, height = 1200)
print(plot_grid(plotlist = p4, rel_widths = c(1, 1.07)))
dev.off()



p3 <- obspred_map(data_cv1_sen, cv1_ml_sen, trans = 'log10', column = FALSE, mask = TRUE)

p4 <- list()

p4[[1]] <- 
  p3[[1]] +
    guides(fill = FALSE) + 
    xlab('Longitude') + 
    ylab('Latitude') +
    ggtitle('Observed data') +
    theme(text = element_text(size = 30)) 



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


png('figs/sen_obs_pred_map_ml.png', width = 1600, height = 1200)
print(plot_grid(plotlist = p4, rel_widths = c(1, 1.07)))
dev.off()




# --------------------------------------------------------------- #
# 2 model comp colombia 
# --------------------------------------------------------------- #



p1 <- obspred_map(data_cv2_col, cv2_covs_col, trans = 'log10',
                  legend_title = 'Cases per 1000',
                  breaks = c(0.01, 0.1, 1, 10, 100), mask = TRUE)
p2 <- obspred_map(data_cv2_col, cv2_ml_col, trans = 'log10', 
                  legend_title = 'API',
                  mask = TRUE)

p1[[1]] <-
  p1[[1]] +
  scale_fill_viridis_c(trans = 'log10', 
                       breaks = c(1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2),
                       labels = trans_format("log10", math_format(10^.x))) 
  

panel1 <- p1[[1]] +
  guides(fill = FALSE) +
  labs(x = '', y = 'Latitude') + 
  lims(x = c(-79.3, -66.5)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))




panel2 <- p1[[2]] +
  guides(fill = FALSE) +
  labs(x = 'Longitude', y = '')+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

panel3 <- p2[[2]] +
  guides(fill = FALSE) +
  labs(x = '', y = '')+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

legend <- get_legend(p1[[1]])

col_preds_plot <- plot_grid(panel1, panel2, panel3, labels = LETTERS[1:3], ncol = 3)
full_plot <- plot_grid(col_preds_plot, legend, ncol = 2, rel_widths = c(7, 1.2))

png('figs/col_comparison_map.png', height = 60, width = 200, unit = 'mm', res = 720)
print(full_plot)
dev.off()



