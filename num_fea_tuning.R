library(dplyr)
library(tidyr)
library(tidyverse)
library(caret)
library(glmnet)
library(e1071)
library(randomForest)
library(pbmcapply)

## training in the discovery cohort and test performance in the validation cohort
load("training_dat.RData") # discovery
load("test_dataset.RData") # validation cohort
load("total_fea.RData") # total features from ##9a. trunk
source("./metabolite_models_tuning.R")

# test 1. non; only metabolite features
# 2.sex, metabolite features with sex
# 3.age, metabolite features with age
# 4.age + sex

set.seed(2024)
num_fea = c(10, 20, 25, 30, 40, 50, 60, 63, 65, 70, 72, 75, 78, 80, 90, 100, 110, 120, 130, 139)
get_results = function(num) {
  n = num_fea[[num]]
  dat = time_fi_all %>%
    mutate(age = `Age.at.assesment..days.`,
           sex_bin = ifelse(Sex == "female", 1, 0)) 
  dat_non = dat %>% 
    dplyr::select(fi, unlist(total_fea[1:n])) %>%
    `colnames<-`(make.names(colnames(.)))
  dat_sex = dat %>% 
    dplyr::select(fi, unlist(total_fea[1:n]), sex_bin) %>% 
    `colnames<-`(make.names(colnames(.)))
  dat_age = dat %>% 
    dplyr::select(fi, unlist(total_fea[1:n]), age) %>% 
    `colnames<-`(make.names(colnames(.)))
  dat_sexage = dat %>% 
    dplyr::select(fi, unlist(total_fea[1:n]), sex_bin, age) %>% 
    `colnames<-`(make.names(colnames(.)))

  rf_model_non = rf_model(dat_non)
  rf_model_sex = rf_model(dat_sex)
  rf_model_age = rf_model(dat_age)
  rf_model_sexage = rf_model(dat_sexage)
  
  nmn_test = nmn_dat %>%
    mutate(age = `Age.at.assesment..days.`,
           sex_bin = ifelse(Sex == "female", 1, 0)) 
  
  nmn_non = nmn_test %>%
    dplyr::select(fi, unlist(total_fea[1:n])) %>% 
    `colnames<-`(make.names(colnames(.)))
  nmn_sex = nmn_test %>% 
    dplyr::select(fi, unlist(total_fea[1:n]), sex_bin) %>% 
    `colnames<-`(make.names(colnames(.)))
  nmn_age = nmn_test %>% 
    dplyr::select(fi, unlist(total_fea[1:n]), age) %>% 
    `colnames<-`(make.names(colnames(.)))
  nmn_sexage = nmn_test %>% 
    dplyr::select(fi, unlist(total_fea[1:n]), sex_bin, age) %>% 
    `colnames<-`(make.names(colnames(.)))
  
  nmn_non_pred = predict(rf_model_non, makeX(nmn_non[, -1]))
  nmn_sex_pred = predict(rf_model_sex, makeX(nmn_sex[, -1]))
  nmn_age_pred = predict(rf_model_age, makeX(nmn_age[, -1]))
  nmn_sexage_pred = predict(rf_model_sexage, makeX(nmn_sexage[, -1]))
  
  rsqured_lst_nmn = sapply(list(unname(nmn_non_pred), 
                                unname(nmn_sex_pred),
                                unname(nmn_age_pred), 
                                unname(nmn_sexage_pred)), function(x) {
                          cor = cor(x, nmn_test$fi)^2
                          return(cor)
                        })
  rmse_lst_nmn = sapply(list(unname(nmn_non_pred), 
                             unname(nmn_sex_pred),
                             unname(nmn_age_pred), 
                             unname(nmn_sexage_pred)), function(x) {
                       rmse = RMSE(x, nmn_test$fi)
                       return(rmse)
                     })
  
  res_df = data.frame(`rsq_nmn` = rsqured_lst_nmn,
                      `rmse_nmn` = rmse_lst_nmn)
  return(res_df)
}

results_102824_2 = pbmclapply(c(1:20), get_results, mc.cores = 20)