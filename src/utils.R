#File containing helper functions, built using updated methodology/tidymodels 
library(tidyverse)
library(tidymodels)
library(caret)

#################
## Data cleaning 

data.generate_samples <- function(df, n){

  #Function that outputs a list of dataframes, with samples of varying size
  #method is meant to mimic scenario where data is collected over time so sample size gets larger.
  #n.vec - vector with samples sizes desired
  #Note function will also return orginal dataframe after removing given rows 
  #Any rows left, are used as the test/holdout set
  
  n.vec <- c(n[1],diff(n))
  list.sample <- list(length(n.vec))
  sample.df <- data.frame(matrix(ncol = ncol(df), nrow = 0))
  for(i in 1:length(n.vec)){
    #get sample and index in original df
    idx.sample <- sample(1:nrow(df),n.vec[i])
    sample.df <- rbind(sample.df,df[idx.sample,])
    #remove those rows from the old dataframe
    df <- df[-idx.sample, ]
    list.sample[[i]] <- sample.df
  }
  
  return(list(list.sample,df))
}

data.clean_data <- function(path){
  #Function that takes a path to a csv and outputs a tidyer dataframe 
  df.raw <- read_csv(path)
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
  
  return(df.full)
}

data.create_matrix <- function(df.full){
  #Function that turns the a dataframe into a model.matrix object 
  #This is because some models require dummy variables for multiple factors 
  
  return(model.matrix( ~ . , df.full))
}

###########################
## Model fitting functions 

#Random forest 
randomforest.fit_model <- function(df, outcome.name= "final_outcome_deathX1"){
  #Function that will fit a random forest with recursive feature elimination 
  #df: training set model.matrix object 
  #outcome.name: outcome in the model 
  #return: fitted model 

  df <- na.omit(df);n <- nrow(df);y.idx <- grep(outcome.name, colnames(df))
  
  #create x and  y matrix 
  df.x <- df[,-y.idx]; df.y <- as.factor(df[,y.idx])
  
  #set control parameters for rfe fitting 
  rfFuncs$summary <- twoClassSummary;
  ctrl <- rfeControl(functions=rfFuncs, method = 'cv', repeats = 5,verbose = F,
                     allowParallel = T,saveDetails = T)
  subsets <- seq(1,25,3)
  rf <- rfe(x=df.x,y=df.y,data = df,sizes = subsets, metric = "ROC", rfeControl = ctrl)
  
  return(rf)
}

randomforest.model_metrics <- function(model,test.df,outcome.name= "final_outcome_deathX1"){
  #Function to compute a range of model metrics, given a model and a test df
  idx <- model$results$Variables == model$optsize;  y.idx <- grep(outcome.name, colnames(test.df))
  roc.val <- model$results[idx,2]
  #compute roc
  r1 <- pROC::roc(test.df[,y.idx], predict(model, test.df, 'prob')[,3]) #computes ROC curve from model applied to test set
  #spec
  ss1 <- data.frame(r1$sensitivities, r1$specificities); names(ss1) <- c('Sens', 'Spec')
  #sens
  spec1 <- dplyr::filter(ss1, abs(Sens - 0.8) == min(abs(Sens - 0.8)))$Spec
  #if more than 1 get mean
  if(length(spec1) > 1){
    spec1 <- mean(spec1,na.rm=T)
  }
  
  #Spec for sens, AUC, Internal ROC
  spec_list <- c(spec1,r1$auc,roc.val)  
  
  #opt variables 
  opt_vars <- model$optVariables
  return(list(spec_list,opt_vars))
}


#Stepwise regression 
stepwise.fit_model <-  function(df, outcome.name= "final_outcome_deathX1"){
  #Function that will fit a stepwise regression model, based on AIC
  #df: training set model.matrix object 
  #outcome.name: outcome in the model 
  #return: fitted model 
  
  df <- na.omit(df);n <- nrow(df);y.idx <- grep(outcome.name, colnames(df))
  
  
  #create x and  y matrix 
  df.x <- df[,-y.idx]; df.y <- as.factor(df[,y.idx])
  #need to set values to yes/no for train function 
  levels(df.y) <- c("no","yes")
  
  train.ctrl <- trainControl("cv", 5,summaryFunction = twoClassSummary,
                             classProbs = T, allowParallel = T,
                             savePredictions = T)
  fit.aic <- train(x = (df.x) , y = df.y,
                   method = "glmStepAIC",family=binomial(),
                   trControl = train.ctrl ,metric = 'ROC', trace=F)
  
  return(fit.aic)
  
}

stepwise.model_metrics <- function(model,test.df,outcome.name= "final_outcome_deathX1"){
  #Function to compute a range of model metrics, given a model and a test df
  y.idx <- grep(outcome.name, colnames(test.df))
  roc.val <- model$results$ROC
  #compute roc
  r1 <- pROC::roc(test.df[,y.idx], predict(model, test.df, 'prob')[,2]) #computes ROC curve from model applied to test set
  #spec
  ss1 <- data.frame(r1$sensitivities, r1$specificities); names(ss1) <- c('Sens', 'Spec')
  #sens
  spec1 <- dplyr::filter(ss1, abs(Sens - 0.8) == min(abs(Sens - 0.8)))$Spec
  #if more than 1 get mean
  if(length(spec1) > 1){
    spec1 <- mean(spec1,na.rm=T)
  }
  
  #Spec for sens, AUC, Internal ROC
  spec_list <- c(spec1,r1$auc,roc.val)  
  #opt variables 
  opt_vars <- variable.names(model$finalModel)
  return(list(spec_list,opt_vars))
}


#Univariate filtering
univariate.fit_model <-  function(df, outcome.name= "final_outcome_deathX1"){
  #Function that will fit a univariate regression model, based on p-value cutoff
  #df: training set model.matrix object 
  #outcome.name: outcome in the model 
  #return: fitted model 
  
  df <- na.omit(df);n <- nrow(df);y.idx <- grep(outcome.name, colnames(df))
  
  
  #create x and  y matrix 
  df.x <- df[,-y.idx]; df.y <- as.factor(df[,y.idx])
  #need to set values to yes/no for train function 
  levels(df.y) <- c("no","yes")
  
  train.ctrl <- trainControl("cv", 5,summaryFunction = twoClassSummary,
                             classProbs = T, allowParallel = T,
                             savePredictions = T)
  fit.glm <- train(x = df.x , y = df.y ,
                   method = "glm",family=binomial(),
                   trControl = train.ctrl,
                   metric = 'ROC', trace=F)
  
  return(fit.glm)
}

univariate.model_metrics <- function(model,test.df,outcome.name= "final_outcome_deathX1"){
  #Function to compute a range of model metrics, given a model and a test df
  y.idx <- grep(outcome.name, colnames(test.df))
  roc.val <- model$results$ROC
  #compute roc
  r1 <- pROC::roc(test.df[,y.idx], predict(model, test.df, 'prob')[,2]) #computes ROC curve from model applied to test set
  #spec
  ss1 <- data.frame(r1$sensitivities, r1$specificities); names(ss1) <- c('Sens', 'Spec')
  #sens
  spec1 <- dplyr::filter(ss1, abs(Sens - 0.8) == min(abs(Sens - 0.8)))$Spec
  #if more than 1 get mean
  if(length(spec1) > 1){
    spec1 <- mean(spec1,na.rm=T)
  }
  
  #Spec for sens ext, AUC ext , Internal ROC
  spec_list <- c(spec1,r1$auc,roc.val)  
  #opt variables 
  coef <- tidy(model$finalModel)
  opt_vars <- coef$term[coef$p.value < 0.05]
  return(list(spec_list,opt_vars))
}


#LASSO 
lasso.fit_model <- function(df, outcome.name= "final_outcome_deathX1"){
  #Function that will fit a LASSO regression model
  #df: training set model.matrix object 
  #outcome.name: outcome in the model 
  #return: fitted model 
  
  df <- na.omit(df);n <- nrow(df);y.idx <- grep(outcome.name, colnames(df))
  
  #create x and  y matrix 
  df.x <- df[,-y.idx]; df.y <- as.factor(df[,y.idx])
  #need to set values to yes/no for train function 
  levels(df.y) <- c("no","yes")
  
  train.ctrl <- trainControl("cv", 5,summaryFunction = twoClassSummary,
                             classProbs = T, allowParallel = T,
                             savePredictions = T)
  fit.glm <- train(x = df.x, y = df.y ,
                   method = 'glmnet', trControl = train.ctrl,
                   metric = 'ROC', trace=F,
                   tuneGrid = expand.grid(alpha = 1,
                                          lambda = seq(0.001,0.1,by = 0.001)),
                   family = "binomial")
  
  return(fit.glm)
}

lasso.model_metrics <- function(model,test.df,outcome.name= "final_outcome_deathX1"){
  #Function to compute a range of model metrics, given a model and a test df
  y.idx <- grep(outcome.name, colnames(test.df))
  roc.val <- max(model$results$ROC)
  #compute roc
  r1 <- pROC::roc(test.df[,y.idx], predict(model, test.df, 'prob')[,2]) #computes ROC curve from model applied to test set
  #spec
  ss1 <- data.frame(r1$sensitivities, r1$specificities); names(ss1) <- c('Sens', 'Spec')
  #sens
  spec1 <- dplyr::filter(ss1, abs(Sens - 0.8) == min(abs(Sens - 0.8)))$Spec
  #if more than 1 get mean
  if(length(spec1) > 1){
    spec1 <- mean(spec1,na.rm=T)
  }
  
  #Spec for sens ext, AUC ext , Internal ROC
  spec_list <- c(spec1,r1$auc,roc.val)  
  final.model <- coef(model$finalModel,model$bestTune$lambda)
  opt_vars <- rownames(final.model)[final.model[,1]!= 0]
  return(list(spec_list,opt_vars))
}


#Elastic net 
en.fit_model <- function(df, outcome.name= "final_outcome_deathX1"){
  #Function that will fit an elastic net regression model
  #df: training set model.matrix object 
  #outcome.name: outcome in the model 
  #return: fitted model 
  
  df <- na.omit(df);n <- nrow(df);y.idx <- grep(outcome.name, colnames(df))
  
  #create x and  y matrix 
  df.x <- df[,-y.idx]; df.y <- as.factor(df[,y.idx])
  #need to set values to yes/no for train function 
  levels(df.y) <- c("no","yes")
  
  tr.ctrl <- trainControl(method="cv", number=5,
                          classProbs=TRUE, summaryFunction=twoClassSummary,
                          allowParallel = T)
  
  fit.glm <- train(x = df.x,y=df.y,
                   method = 'glmnet', trControl = tr.ctrl,
                   metric = 'ROC', trace=F,
                   tuneGrid = expand.grid(alpha = seq(0.1, 1, length = 10),
                                          lambda = seq(0.001,0.1,by = 0.001)),
                   family = "binomial")
  return(fit.glm)
}

en.model_metrics <- function(model,test.df,outcome.name= "final_outcome_deathX1"){
  #Function to compute a range of model metrics, given a model and a test df
  y.idx <- grep(outcome.name, colnames(test.df))
  roc.val <- max(model$results$ROC)
  #compute roc
  r1 <- pROC::roc(test.df[,y.idx], predict(model, test.df, 'prob')[,2]) #computes ROC curve from model applied to test set
  #spec
  ss1 <- data.frame(r1$sensitivities, r1$specificities); names(ss1) <- c('Sens', 'Spec')
  #sens
  spec1 <- dplyr::filter(ss1, abs(Sens - 0.8) == min(abs(Sens - 0.8)))$Spec
  #if more than 1 get mean
  if(length(spec1) > 1){
    spec1 <- mean(spec1,na.rm=T)
  }
  
  #Spec for sens ext, AUC ext , Internal ROC
  spec_list <- c(spec1,r1$auc,roc.val)  
  final.model <- coef(model$finalModel,model$bestTune$lambda)
  opt_vars <- rownames(final.model)[final.model[,1]!= 0]
  return(list(spec_list,opt_vars))
}


#Adapative LASSO
adaLasso.fit_model <- function(df, outcome.name= "final_outcome_deathX1"){
  #Function that will fit a LASSO regression model
  #df: training set model.matrix object 
  #outcome.name: outcome in the model 
  #return: fitted model 
  
  df <- na.omit(df);n <- nrow(df);y.idx <- grep(outcome.name, colnames(df))
  
  #create x and  y matrix 
  df.x <- df[,-y.idx]; df.y <- as.factor(df[,y.idx])
  #need to set values to yes/no for train function 
  levels(df.y) <- c("no","yes")
  
  tr.ctrl <- trainControl(method="cv", number=5,
                          classProbs=TRUE, summaryFunction=twoClassSummary,
                          allowParallel = T)
  #fit ridge regression
  fit.glm <- train(x = df.x,y=df.y,
                   method = 'glmnet', trControl = tr.ctrl,
                   metric = 'ROC', trace=F,
                   tuneGrid = expand.grid(alpha = 0,
                                          lambda = seq(0.001,0.1,by = 0.001)),
                   family = "binomial")
  
  #get coefs of the vars @ best lambda
  best.ridge.coefs <- as.matrix(coef(fit.glm$finalModel,fit.glm$bestTune$lambda))[-1]
  
  #adalassao
  fit.adaLasso <- train(x = df.x,y=df.y,
                        method = 'glmnet', trControl = tr.ctrl,
                        metric = 'ROC', trace=F,
                        tuneGrid = expand.grid(alpha = 1,
                                               lambda = seq(0.001,0.1,by = 0.001)),
                        family = "binomial" ,
                        penalty.factor = 1/abs(best.ridge.coefs))
  
  
  return(fit.glm)
}

adaLasso.model_metrics <- function(model,test.df,outcome.name= "final_outcome_deathX1"){
  #Function to compute a range of model metrics, given a model and a test df
  y.idx <- grep(outcome.name, colnames(test.df))
  roc.val <- max(model$results$ROC)
  #compute roc
  r1 <- pROC::roc(test.df[,y.idx], predict(model, test.df, 'prob')[,2]) #computes ROC curve from model applied to test set
  #spec
  ss1 <- data.frame(r1$sensitivities, r1$specificities); names(ss1) <- c('Sens', 'Spec')
  #sens
  spec1 <- dplyr::filter(ss1, abs(Sens - 0.8) == min(abs(Sens - 0.8)))$Spec
  #if more than 1 get mean
  if(length(spec1) > 1){
    spec1 <- mean(spec1,na.rm=T)
  }
  
  #Spec for sens ext, AUC ext , Internal ROC
  spec_list <- c(spec1,r1$auc,roc.val)  
  final.model <- coef(model$finalModel,model$bestTune$lambda)
  opt_vars <- rownames(final.model)[final.model[,1]!= 0]
  return(list(spec_list,opt_vars))
}

#RFE 
rfe.fit_model <- function(df, outcome.name= "final_outcome_deathX1"){
  #Function that will fit a random forest with recursive feature elimination 
  #df: training set model.matrix object 
  #outcome.name: outcome in the model 
  #return: fitted model 
  
  df <- na.omit(df);n <- nrow(df);y.idx <- grep(outcome.name, colnames(df))
  
  #create x and  y matrix 
  df.x <- df[,-y.idx]; df.y <- as.factor(df[,y.idx])
  
  #set control parameters for rfe fitting 
  lrFuncs$summary <- twoClassSummary;
  ctrl <- rfeControl(functions=lrFuncs, method = 'cv', repeats = 5,verbose = F, allowParallel = T,saveDetails = T)
  subsets <- seq(1,25,3)
  rfe.model <- rfe(x=df.x,y=df.y,sizes = subsets, metric = "ROC", rfeControl = ctrl)
  
  return(rfe.model)
}

rfe.model_metrics <- function(model,test.df,outcome.name= "final_outcome_deathX1"){
  #Function to compute a range of model metrics, given a model and a test df
  idx <- model$results$Variables == model$optsize;  y.idx <- grep(outcome.name, colnames(test.df))
  roc.val <- model$results[idx,2]
  #compute roc
  r1 <- pROC::roc(test.df[,y.idx], predict(model, test.df, 'prob')[,2]) #computes ROC curve from model applied to test set
  #spec
  ss1 <- data.frame(r1$sensitivities, r1$specificities); names(ss1) <- c('Sens', 'Spec')
  #sens
  spec1 <- dplyr::filter(ss1, abs(Sens - 0.8) == min(abs(Sens - 0.8)))$Spec
  #if more than 1 get mean
  if(length(spec1) > 1){
    spec1 <- mean(spec1,na.rm=T)
  }
  
  #Spec for sens, AUC, Internal ROC
  spec_list <- c(spec1,r1$auc,roc.val)  
  #opt variables 
  opt_vars <- model$optVariables
  return(list(spec_list,opt_vars))
}


#######
#Misc

add_to_matrix <- function(model.list,metrics.list,col_names,n,frequency.matrix,method){
  #Method ot add row to frequency matrix
  
  for(i in 1:length(metrics.list)){
    top_vars <- metrics.list[[i]][[2]]
    in.model <- as.numeric(col_names %in% top_vars)
    model.metric <- metrics.list[[i]][[1]]
    row <- c(in.model,method,n[i],model.metric)
    frequency.matrix[nrow(frequency.matrix) + 1,] = row
  }
  
  return(frequency.matrix)
}







