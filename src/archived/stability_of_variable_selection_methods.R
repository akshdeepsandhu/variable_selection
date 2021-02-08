#load packages
library(tidyverse)
library(glmnet)
library(tableone)
library(caret)
library(ROCR)
library(parallel)
library(broom)
library(tictoc)

cat('Starting...\n')
#load and reformat data
df.raw <- read_csv('../data/SD_for_Ash.csv')
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
#model fitting procedure


do_model_fitting <- function(df.full,df.matrix, rep.val){

  #(a): set up matricies for data collection

  model.vars <- colnames(df.matrix)[-1]
  n <- c(500,1000,2000,2500,3000,3500,4000)
  method <- c("RF", "stepwiseAIC", "univarate.pval", "LASSO", "RFE", "EN")

  #variable selection frequency matrix
  frequency.matrix <- data.frame(matrix(ncol = ncol(df.matrix) + 2, nrow = 0))
  colnames(frequency.matrix) <- c(model.vars,"method", "n", "ROC")

  #helper function to generate samples of size n
  get_sample <- function(df, n){
    #outputs a list of dataframes, with samples of varying size
    #method is meant to mimic scenario where data is collected over time so sample size gets larger.
    #n.vec - vector with samples sizes desired

    n.vec <- c(n[1],diff(n))
    list.sample <- list(length(n.vec))
    sample.df <- data.frame(matrix(ncol = ncol(df), nrow = 0))
    for(i in 1:length(n.vec)){
      #get sample and index in original df
      sample.df <- rbind(sample.df,sample_n(as.data.frame(df),size = n.vec[i]))
      idx.sample <- plyr::match_df(as.data.frame(df),sample.df) %>% rownames() %>% as.numeric()
      #remove those rows from the old dataframe
      df <- df[-idx.sample, ]
      list.sample[[i]] <- sample.df
    }

    return(list.sample)
  }

  #(b): model fitting methods

  #rf method
  get_rf_results <- function(sample.list, col.names){



    #helper to speed up
    get_rf_vars <- function(sample.df, names.vec){
      #Get the top variables chosen by RF for different samples sizes

      #helper functions
      fit.rf <- function(df){
        #function that will fit a rf and return top 10 most important varibales

        #preprocess data for rf fitting
        df.matrix <- model.matrix(final_outcome_death ~ . , df)
        n <- nrow(df); model.vars <- colnames(df.matrix)[-1]
        df <- na.omit(df)
        x.mat <- df.matrix[,-1]
        df$final_outcome_death <- factor(make.names(df$final_outcome_death))
        #fit rf
        rfFuncs$summary <- twoClassSummary;
        ctrl <- rfeControl(functions=rfFuncs, method = 'cv', repeats = 5,verbose = F,
                           allowParallel = T,saveDetails = T)
        subsets <- seq(1,25,3)
        rf <- rfe(x = x.mat, y = df$final_outcome_death, sizes = subsets, metric = "ROC", rfeControl = ctrl)
        #get top variables
        in.model <- as.numeric(model.vars %in% rf$optVariables)
        #get roc value
        idx <- rf$results$Variables == rf$optsize
        roc.val <- rf$results[idx,2]


        return(c(in.model,"RF",n, roc.val))

      }

      #fit row and get vars chosen
      out.row <- as.data.frame(t(fit.rf(sample.df)))
      names(out.row) <- names.vec
      return(out.row)

    }

    #fit in parallel and get the results
    results <- mclapply(sample.list, get_rf_vars,
                        names.vec = col.names, mc.cores = detectCores())
    #Fitting takes ~ 480seconds

    return(results)
  }

  #stepwise
  get_stepwise_results <- function(sample.list, col.names){



    #helper to speed up
    get_stepaic_vars <- function(sample.df, names.vec){

      fit.stepwise <- function(df){

        #create model matrix for fitting
        df.matrix <- model.matrix(final_outcome_death ~ . , df)
        #get n and variable names
        n <- nrow(df); model.vars <- colnames(df.matrix)[-1]
        #clean up dataframe to allow for model fitting
        df <- df %>% drop_na()
        #fit stepwiseAIC regression to glm object
        x.mat <- df.matrix[,-1]
        train.ctrl <- trainControl("cv", 5,summaryFunction = twoClassSummary,
                                   classProbs = T, allowParallel = T,
                                   savePredictions = T)
        fit.aic <- train(x = x.mat , y = df$final_outcome_death,
                         method = "glmStepAIC",family=binomial(),
                         trControl = train.ctrl ,metric = 'ROC', trace=F)
        #get vars
        top.vars <- variable.names(fit.aic$finalModel)[-1]
        #score variable with 1 if it is in the model and 0 o.w
        in.model <- as.numeric(model.vars %in% top.vars)
        roc.val <- fit.aic$results$ROC
        return(c(in.model,"stepAIC",n, roc.val))

      }

      #fit row and get vars chosen
      out.row <- as.data.frame(t(fit.stepwise(sample.df)))
      names(out.row) <- names.vec
      return(out.row)

    }


    #fit in parallel and get the results
    results <- mclapply(sample.list, get_stepaic_vars,
                        names.vec = col.names,
                        mc.cores = detectCores())

    return(results)

  }

  #univariate
  get_univarate_results <- function(sample.list, col.names){




    #helper to speed up
    get_univariate_vars <- function(sample.df, names.vec){

      fit.univarate <- function(df){

        #create model matrix for fitting
        df.matrix <- model.matrix(final_outcome_death ~ . , df)
        #get n and variable names
        n <- nrow(df); model.vars <- colnames(df.matrix)[-1]
        #clean up dataframe to allow for model fitting
        df <- df %>% drop_na()
        #fit stepwiseAIC regression to glm object
        x.mat <- df.matrix[,-1]
        ctrl <- trainControl("cv", 5,classProbs=TRUE, summaryFunction = twoClassSummary,
                             allowParallel = T)
        fit.glm <- train(x = x.mat , y = df$final_outcome_death,
                         method = "glm",family=binomial(),
                         trControl = ctrl,
                         metric = 'ROC', trace=F)
        #get vars
        coef <- tidy(fit.glm$finalModel)
        coef <- coef[-1,]
        top.vars <- coef$term[coef$p.value < 0.05]
        #score variable with 1 if it is in the model and 0 o.w
        in.model <- as.numeric(model.vars %in% top.vars)
        roc.val <- fit.glm$results$ROC
        return(c(in.model,"univariate.pval",n, roc.val))

      }

      #fit row and get vars chosen
      out.row <- as.data.frame(t(fit.univarate(sample.df)))
      names(out.row) <- names.vec
      return(out.row)

    }


    #fit in parallel and get the results
    results <- mclapply(sample.list, get_univariate_vars,
                        names.vec = col.names,
                        mc.cores = detectCores())

    return(results)

  }

  #lasso
  get_lasso_results <- function(sample.list, col.names){




    #helper to speed up
    get_lasso_vars <- function(sample.df, names.vec){

      fit.lasso <- function(df){

        #create model matrix for fitting
        df.matrix <- model.matrix(final_outcome_death ~ . , df)
        #get n and variable names
        n <- nrow(df); model.vars <- colnames(df.matrix)[-1]
        #clean up dataframe to allow for model fitting
        df <- df %>% drop_na()
        #fit stepwiseAIC regression to glm object
        x.mat <- df.matrix[,-1]
        #fit lasso
        tr.ctrl <- trainControl(method="cv", number=5,
                     classProbs=TRUE, summaryFunction=twoClassSummary,
                     allowParallel = T)

        fit.glm <- train(x = x.mat,y=df$final_outcome_death,
                         method = 'glmnet', trControl = tr.ctrl,
                         metric = 'ROC', trace=F,
                         tuneGrid = expand.grid(alpha = 1,
                                                lambda = seq(0.001,0.1,by = 0.001)),
                         family = "binomial")

        #get vars
        final.model <- coef(fit.glm$finalModel,fit.glm$bestTune$lambda)
        top.vars <- rownames(final.model)[final.model[,1]!= 0]
        #score variable with 1 if it is in the model and 0 o.w
        roc.val <- max(fit.glm$results$ROC)
        in.model <- as.numeric(model.vars %in% top.vars)

        return(c(in.model,"LASSO",n, roc.val))

      }

      #fit row and get vars chosen
      out.row <- as.data.frame(t(fit.lasso(sample.df)))
      names(out.row) <- names.vec
      return(out.row)

    }


    #fit in parallel and get the results
    results <- mclapply(sample.list, get_lasso_vars,
                        names.vec = col.names,
                        mc.cores = detectCores())

    return(results)

  }

  #rfe with logistic regression
  get_rfe_results <- function(sample.list, col.names){



    #helper to speed up
    get_rfe_vars <- function(sample.df, names.vec){
      #Get the top variables chosen by RFE for different samples sizes

      #helper functions
      fit.rfe <- function(df){
        #function that will fit a rfe and return top 10 most important varibales

        #preprocess data for rf fitting
        df.matrix <- model.matrix(final_outcome_death ~ . , df)
        n <- nrow(df); model.vars <- colnames(df.matrix)[-1]
        df <- na.omit(df)
        x.mat <- df.matrix[,-1]
        df$final_outcome_death <- factor(make.names(df$final_outcome_death))
        #fit rf
        lrFuncs$summary <- twoClassSummary;
        ctrl <- rfeControl(functions=lrFuncs, method = 'cv', repeats = 5,verbose = F, allowParallel = T,saveDetails = T)
        subsets <- seq(1,25,3)
        rfe.model <- rfe(x = x.mat, y = df$final_outcome_death, sizes = subsets, metric = "ROC", rfeControl = ctrl)
        #get top variables
        in.model <- as.numeric(model.vars %in% rfe.model$optVariables)
        #get roc value
        idx <- rfe.model$results$Variables == rfe.model$optsize
        roc.val <- rfe.model$results[idx,2]


        return(c(in.model,"RFE",n, roc.val))

      }

      #fit row and get vars chosen
      out.row <- as.data.frame(t(fit.rfe(sample.df)))
      names(out.row) <- names.vec
      return(out.row)

    }

    #fit in parallel and get the results
    results <- mclapply(sample.list, get_rfe_vars,
                        names.vec = col.names, mc.cores = detectCores())
    #Fitting takes ~ 480seconds

    return(results)
  }

  #elastic net
  get_en_results <- function(sample.list, col.names){




    #helper to speed up
    get_en_vars <- function(sample.df, names.vec){

      fit.en <- function(df){

        #create model matrix for fitting
        df.matrix <- model.matrix(final_outcome_death ~ . , df)
        #get n and variable names
        n <- nrow(df); model.vars <- colnames(df.matrix)[-1]
        #clean up dataframe to allow for model fitting
        df <- df %>% drop_na()
        #fit stepwiseAIC regression to glm object
        x.mat <- df.matrix[,-1]
        #fit lasso
        tr.ctrl <- trainControl(method="cv", number=5,
                                classProbs=TRUE, summaryFunction=twoClassSummary,
                                allowParallel = T)

        fit.glm <- train(x = x.mat,y=df$final_outcome_death,
                         method = 'glmnet', trControl = tr.ctrl,
                         metric = 'ROC', trace=F,
                         tuneGrid = expand.grid(alpha = seq(0.1, 1, length = 10),
                                                lambda = seq(0.001,0.1,by = 0.001)),
                         family = "binomial")

        #get vars
        final.model <- coef(fit.glm$finalModel,fit.glm$bestTune$lambda)
        top.vars <- rownames(final.model)[final.model[,1]!= 0]
        #score variable with 1 if it is in the model and 0 o.w
        roc.val <- max(fit.glm$results$ROC)
        in.model <- as.numeric(model.vars %in% top.vars)

        return(c(in.model,"EN",n, roc.val))

      }

      #fit row and get vars chosen
      out.row <- as.data.frame(t(fit.en(sample.df)))
      names(out.row) <- names.vec
      return(out.row)

    }


    #fit in parallel and get the results
    results <- mclapply(sample.list, get_en_vars,
                        names.vec = col.names,
                        mc.cores = detectCores())

    return(results)

  }

  #(c): fit models to samples

  #generate samples
  sample.list <- get_sample((df.full),n)
  columns <- colnames(frequency.matrix)
  #get rf results
  rf.results <- get_rf_results(sample.list, columns)
  frequency.matrix <- do.call(rbind,rf.results)
  cat("Finished getting RF results\n")
  #get stepwise results
  stepwise.results <- get_stepwise_results(sample.list,columns)
  stepwise.results <- do.call(rbind,stepwise.results)
  frequency.matrix <- rbind(frequency.matrix,stepwise.results)
  ("Finished getting stepwiseAIC results\n")
  #get univariate resutls
  univariate.results <- get_univarate_results(sample.list,columns)
  univariate.results <- do.call(rbind,univariate.results)
  frequency.matrix <- rbind(frequency.matrix, univariate.results)
  cat("Finished getting univariate filtering results\n")
  #get LASSO resutls
  lasso.results <- get_lasso_results(sample.list,columns)
  lasso.results <- do.call(rbind,lasso.results)
  frequency.matrix <- rbind(frequency.matrix, lasso.results)
  cat("Finished getting LASSO results\n")
  #rfe with lr
  rfeLr.results <- get_rfe_results(sample.list,columns)
  rfeLr.results <- do.call(rbind,rfeLr.results)
  frequency.matrix <- rbind(frequency.matrix,rfeLr.results)
  cat("Finished getting RFE results\n")
  #get en results
  en.results <- get_en_results(sample.list,columns)
  en.results <- do.call(rbind,en.results)
  frequency.matrix <- rbind(frequency.matrix, en.results)
  cat("Finished getting EN results\n")
  #add iteration
  frequency.matrix$rep_val <- rep.val
  cat("----------------------------\n")
  return(frequency.matrix)
}


#1 run takes ~ 620 seconds
out <- list(500)
for(i in 1:500){

  out[[i]] <- do_model_fitting(df.full,df.matrix,i)
  cat('Completed' ,i,  'bootsrap run\n' )
}

cat('Writing results to file\n')
out <- do.call(rbind,out)
write.csv(out,'~/Desktop/results_21_2_20_500_runs.csv', row.names = F)
cat('Done')
