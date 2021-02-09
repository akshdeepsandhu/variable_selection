#load packages
library(tidyverse)
library(glmnet)
library(tableone)
library(caret)
library(ROCR)
library(parallel)
library(broom)
library(tictoc)
library(lme4)
library(lmerTest)
source('src/fit_funcs.R')
cat('Starting...\n')

#BELOW SCRIPT WORKS: 28/03/2020
################################
#process data
#load and reformat data
#path  = 'Smart Discharge/data/SD_for_Ash.csv'
df.raw <- read_csv('Smart Discharge/data/SD_for_Ash.csv')
#change type to categorical for relevant variables
factor.vars <- c("sex","Abnormal_BCS","hivFinal","bednet","final_outcome_death", "boiling", "hivMaternal", "malaria", "sibling_deaths")
del.vars <- c("water_source_drinking")
#clean up
df.raw <- df.raw %>%
  select(-del.vars) %>%
  mutate_at(vars(factor.vars), funs(factor))
#change factor to str
df.raw$final_outcome_death <- make.names(df.raw$final_outcome_death)
#centre and scale
df.scaled <- scale(df.raw %>% select(-c(factor.vars, W_AZ, W_LZ,BMI_AZ)))
df.full <- cbind(as.data.frame(df.scaled), df.raw %>% select(factor.vars, W_AZ, W_LZ, BMI_AZ))
df.full <- df.full %>%
  mutate(temp2 = temperature*temperature)
df.matrix <- model.matrix(final_outcome_death ~ . , df.full)
#remove df
rm(df.raw)
###################################

do_model_fitting <- function(df.full,df.matrix, rep.val){

  #(a): set up matricies for data collection

  model.vars <- colnames(df.matrix)[-1]
  n <- c(500,1000,2000,2500,3000,3500,4000)
  method <- c("RF", "stepwiseAIC", "univarate.pval", "LASSO", "RFE", "EN")

  #variable selection frequency matrix
  frequency.matrix <- data.frame(matrix(ncol = ncol(df.matrix) + 2, nrow = 0))
  colnames(frequency.matrix) <- c(model.vars,"method", "N", "ROC")

  #(b): fit models to samples

  #generate samples
  sample.out <- get_sample((df.full),n)
  sample.list <- sample.out[[1]]
  test.df <- sample.out[[2]]
  column.names <- colnames(frequency.matrix)[-c(29:31)]
  y.name <- "final_outcome_death"

  #get rf results
  rf.results <- get_rf_results(sample.list, column.names,y.name,test.df)
  frequency.matrix <- do.call(rbind,rf.results)
  cat("Finished getting RF results\n")
  #get stepwise results
  stepwise.results <- get_stepwise_results(sample.list, column.names,y.name,test.df)
  stepwise.results <- do.call(rbind,stepwise.results)
  frequency.matrix <- rbind(frequency.matrix,stepwise.results)
  ("Finished getting stepwiseAIC results\n")
  #get univariate resutls
  univariate.results <- get_univarate_results(sample.list,column.names,y.name,test.df)
  univariate.results <- do.call(rbind,univariate.results)
  frequency.matrix <- rbind(frequency.matrix, univariate.results)
  cat("Finished getting univariate filtering results\n")
  #get LASSO resutls
  lasso.results <- get_lasso_results(sample.list,column.names,y.name,test.df)
  lasso.results <- do.call(rbind,lasso.results)
  frequency.matrix <- rbind(frequency.matrix, lasso.results)
  cat("Finished getting LASSO results\n")
  #rfe with lr
  rfeLr.results <- get_rfe_results(sample.list,column.names,y.name,test.df)
  rfeLr.results <- do.call(rbind,rfeLr.results)
  frequency.matrix <- rbind(frequency.matrix,rfeLr.results)
  cat("Finished getting RFE results\n")
  #get en results
  en.results <- get_en_results(sample.list,column.names,y.name,test.df)
  en.results <- do.call(rbind,en.results)
  frequency.matrix <- rbind(frequency.matrix, en.results)
  cat("Finished getting EN results\n")
  #get adaLasso results
  adaLasso.results <- get_adaLasso_results(sample.list,column.names,y.name,test.df)
  adaLasso.results <- do.call(rbind,adaLasso.results)
  frequency.matrix <- rbind(frequency.matrix,adaLasso.results)

  #add iteration
  frequency.matrix$rep_val <- rep.val
  cat("----------------------------\n")
  return(frequency.matrix)
  }

################################################
#1 run takes ~ 620 seconds
out <- list(500)
for(i in 1:500){

  out[[i]] <- do_model_fitting(df.full,df.matrix,i)

  cat('Completed' ,i,  'bootsrap run\n' )
}


cat('Writing results to file\n')
out <- do.call(rbind,out)
write.csv(x=as.data.frame(out), file = 'resultst_test_run.csv', row.names = F)
cat('Done')


####New work
df.matrix <- data.create_matrix(data.clean_data('data/SD_for_Ash.csv'))



















