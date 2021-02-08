#FittingFuncs.R
#R file containing functions that will be used for model fitting


#sample generation method
get_sample <- function(df, n){
  #outputs a list of dataframes, with samples of varying size
  #method is meant to mimic scenario where data is collected over time so sample size gets larger.
  #n.vec - vector with samples sizes desired

  n.vec <- c(n[1],diff(n))
  list.sample <- list(length(n.vec))
  sample.df <- data.frame(matrix(ncol = ncol(df), nrow = 0))
  print(n.vec)
  for(i in 1:length(n.vec)){
    #get sample and index in original df
    idx.sample <- sample(1:nrow(df),n.vec[i])
    sample.df <- rbind(sample.df,df[idx.sample,])
    print(nrow(sample.df))
    #remove those rows from the old dataframe
    df <- df[-idx.sample, ]
    print(nrow(df))
    list.sample[[i]] <- sample.df
  }

  return(list(list.sample,df))
}


#stepwise
get_rf_results <- function(sample.list, model.vars, outcome.name, test.df){
  #Fits in parallel random forest rfe to a list of dataframes
  #param: sample.list: list of dataframes containing sample data
  #param: model.vars: the variable in dataframe
  #param: outcome.name: string name of outcome

  spec_for_sens <- function(model,test.df){
    #compute roc
    r1 <- pROC::roc(test.df$final_outcome_death, predict(model, test.df, 'prob')$X1) #computes ROC curve from model applied to test set
    #spec
    ss1 <- data.frame(r1$sensitivities, r1$specificities); names(ss1) <- c('Sens', "Spec")
    #sens
    spec1 <- dplyr::filter(ss1, abs(Sens - 0.8) == min(abs(Sens - 0.8)))$Spec
    #if more than 1 get mean
    if(length(spec1) > 1){
      spec1 <- mean(spec1,na.rm=T)
    }
    return(spec1)

  }


  #helper to speed up
  get_rf_vars <- function(sample.df, model.vars,outcome.name,test.df){
    #Get the top variables chosen by RF for different samples sizes

    #helper functions
    fit.rf <- function(df, model.vars,outcome.name,test.df){
      #function that will fit a rf and return top 10 most important varibales
      #get n and variable names
      n <- nrow(df); df <- na.omit(df); y.idx <- grep(outcome.name, colnames(df))
      #get x and y dataframes
      df.x <- df[,-y.idx]; df.y <- as.factor(df[,y.idx])
      #fit rf
      rfFuncs$summary <- twoClassSummary;
      ctrl <- rfeControl(functions=rfFuncs, method = 'cv', repeats = 5,verbose = F,
                         allowParallel = T,saveDetails = T)
      subsets <- seq(1,25,3)
      rf <- rfe(x=df.x,y=df.y,data = df,sizes = subsets, metric = "ROC", rfeControl = ctrl)
      #get top variables
      in.model <- as.numeric(model.vars %in% rf$optVariables)
      #get roc value
      idx <- rf$results$Variables == rf$optsize
      roc.val <- rf$results[idx,2]

      spec.sens <- spec_for_sens(rf,test.df)

      return(c(in.model,"RF",n, roc.val, spec.sens))

    }

    #fit row and get vars chosen
    out.row <- as.data.frame(t(fit.rf(sample.df, model.vars,outcome.name,test.df)))
    model.vars <- c(model.vars, "method", "N", "ROC", "spec_for_sens")
    names(out.row) <- model.vars
    return(out.row)

  }

  #fit in parallel and get the results
  #make cluster
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, {library(caret)})
  results <- parLapply(cl,sample.list, get_rf_vars,
                       model.vars = model.vars, outcome.name = outcome.name,
                       test.df=test.df)
  stopCluster(cl)
  #Fitting takes ~ 480seconds

  return(results)
}

#stepwise
get_stepwise_results <- function(sample.list, model.vars,outcome.name, test.df){
  #Fits in parallel stepwiseAIC to a list of dataframes
  #param: sample.list: list of dataframes containing sample data
  #param: model.vars: the variable in dataframe
  #param: outcome.name: string name of outcome

  spec_for_sens <- function(model,test.df){
    #compute roc
    r1 <- pROC::roc(test.df$final_outcome_death, predict(model, test.df, 'prob')$X1) #computes ROC curve from model applied to test set
    #spec
    ss1 <- data.frame(r1$sensitivities, r1$specificities); names(ss1) <- c('Sens', "Spec")
    #sens
    spec1 <- dplyr::filter(ss1, abs(Sens - 0.8) == min(abs(Sens - 0.8)))$Spec
    #if more than 1 get mean
    if(length(spec1) > 1){
      spec1 <- mean(spec1,na.rm=T)
    }
    return(spec1)

  }

  #helper to speed up
  get_stepaic_vars <- function(sample.df, model.vars,outcome.name,test.df){

    fit.stepwise <- function(df, model.vars,outcome.name,test.df){

      #get n and variable names
      n <- nrow(df); df <- na.omit(df); y.idx <- grep(outcome.name, colnames(df))
      #get x and y dataframes
      df.x <- df[,-y.idx]; df.y <- as.factor(df[,y.idx])
      #fit stepwiseAIC regression to glm object
      train.ctrl <- trainControl("cv", 5,summaryFunction = twoClassSummary,
                                 classProbs = T, allowParallel = T,
                                 savePredictions = T)
      fit.aic <- train(x = df.x , y = df.y,
                       method = "glmStepAIC",family=binomial(),
                       trControl = train.ctrl ,metric = 'ROC', trace=F)

      #get vars
      top.vars <- variable.names(fit.aic$finalModel)
      #score variable with 1 if it is in the model and 0 o.w
      in.model <- as.numeric(model.vars %in% top.vars)
      roc.val <- fit.aic$results$ROC

      spec.sens <- spec_for_sens(fit.aic,test.df)

      return(c(in.model,"stepAIC",n, roc.val, spec.sens))

    }

    #fit row and get vars chosen
    out.row <- as.data.frame(t(fit.stepwise(sample.df,model.vars,outcome.name,test.df)))
    model.vars <-c(model.vars, "method", "N", "ROC", "spec_for_sens")
    names(out.row) <- model.vars
    return(out.row)


  }


  #fit in parallel and get the results
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, library(caret))
  results <- parLapply(cl, sample.list, get_stepaic_vars,
                       model.vars = model.vars, outcome.name = outcome.name,
                       test.df = test.df)
  stopCluster(cl)

  return(results)

}

#univariate
get_univarate_results <- function(sample.list, model.vars,outcome.name,test.df){
  #Fits in parallel univariate p-value selection to a list of dataframes
  #param: sample.list: list of dataframes containing sample data
  #param: model.vars: the variable in dataframe
  #param: outcome.name: string name of outcome

  spec_for_sens <- function(model,test.df){
    #compute roc
    r1 <- pROC::roc(test.df$final_outcome_death, predict(model, test.df, 'prob')$X1) #computes ROC curve from model applied to test set
    #spec
    ss1 <- data.frame(r1$sensitivities, r1$specificities); names(ss1) <- c('Sens', "Spec")
    #sens
    spec1 <- dplyr::filter(ss1, abs(Sens - 0.8) == min(abs(Sens - 0.8)))$Spec
    #if more than 1 get mean
    if(length(spec1) > 1){
      spec1 <- mean(spec1,na.rm=T)
    }
    return(spec1)

  }
  #helper to speed up
  get_univariate_vars <- function(sample.df, model.vars,outcome.name,test.df){

    fit.univarate <- function(df,model.vars,outcome.name,test.df){

      #get n and variable names
      n <- nrow(df); df <- na.omit(df); y.idx <- grep(outcome.name,colnames(df))
      #get x and y dataframes
      df.x <- df[,-y.idx]; df.y <- as.factor(df[,y.idx])
      #fit lm
      train.ctrl <- trainControl("cv", 5,summaryFunction = twoClassSummary,
                                 classProbs = T, allowParallel = T,
                                 savePredictions = T)
      fit.glm <- train(x = df.x , y = df.y ,
                       method = "glm",family=binomial(),
                       trControl = train.ctrl,
                       metric = 'ROC', trace=F)
      #get vars
      coef <- tidy(fit.glm$finalModel)
      coef <- coef[-1,]
      top.vars <- coef$term[coef$p.value < 0.05]
      #score variable with 1 if it is in the model and 0 o.w
      in.model <- as.numeric(model.vars %in% top.vars)
      roc.val <- fit.glm$results$ROC

      spec.sens <- spec_for_sens(fit.glm,test.df)

      return(c(in.model,"univariate.pval",n, roc.val,spec.sens))

    }

    #fit row and get vars chosen
    out.row <- as.data.frame(t(fit.univarate(sample.df,model.vars,outcome.name,test.df)))
    model.vars <- c(model.vars, "method", "N", "ROC", "spec_for_sens")
    names(out.row) <- model.vars
    return(out.row)

  }


  #fit in parallel and get the results
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, {library(caret)
    library(broom)})
  results <- parLapply(cl, sample.list, get_univariate_vars,
                       model.vars = model.vars,outcome.name=outcome.name,
                       test.df =  test.df)
  stopCluster(cl)
  return(results)


}

#lasso
get_lasso_results <- function(sample.list, model.vars,outcome.name,test.df){
  #Fits in parallel LASSO regressiom for variable selectio to a list of dataframes
  #param: sample.list: list of dataframes containing sample data
  #param: model.vars: the variable in dataframe
  #param: outcome.name: string name of outcome
  spec_for_sens <- function(model,outcome.name,test.df){
    newx = model.matrix(formula(paste(outcome.name, "~", ".")), data = test.df)
    #compute roc
    lasso_prob <- predict(model, type="prob",
                          newdata  = newx , s = model$finalModel$lambdaOpt)

    r1 <- pROC::roc(test.df$final_outcome_death, lasso_prob$X1)  #computes ROC curve from model applied to test set
    #spec
    ss1 <- data.frame(r1$sensitivities, r1$specificities); names(ss1) <- c('Sens', "Spec")
    #sens
    spec1 <- dplyr::filter(ss1, abs(Sens - 0.8) == min(abs(Sens - 0.8)))$Spec
    #if more than 1 get mean
    if(length(spec1) > 1){
      spec1 <- mean(spec1,na.rm=T)
    }
    return(spec1)

  }
  #helper to speed up
  get_lasso_vars <- function(sample.df, model.vars,outcome.name,test.df){

    fit.lasso <- function(df, model.vars,outcome.name,test.df){


      #get n and variable names
      n <- nrow(df); df <- na.omit(df); y.idx <- grep(outcome.name,colnames(df))
      #get x and y dataframes
      df.x <- model.matrix(formula(paste(outcome.name, "~", ".")), data = df); df.y <- as.factor(df[,y.idx])
      #fit LASSO
      train.ctrl <- trainControl("cv", 5,summaryFunction = twoClassSummary,
                                 classProbs = T, allowParallel = T,
                                 savePredictions = T)
      fit.glm <- train(x = df.x, y = df.y ,
                       method = 'glmnet', trControl = train.ctrl,
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

      spec.sens <- spec_for_sens(fit.glm,outcome.name,test.df)

      return(c(in.model,"LASSO",n, roc.val, spec.sens))

    }

    #fit row and get vars chosen
    out.row <- as.data.frame(t(fit.lasso(sample.df,model.vars,outcome.name,test.df)))
    model.vars <- c(model.vars,  "method", "N", "ROC", "spec_for_sens")
    names(out.row) <- model.vars
    return(out.row)

  }


  #fit in parallel and get the results
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, library(caret))
  results <- parLapply(cl, sample.list, get_lasso_vars,
                       model.vars = model.vars, outcome.name=outcome.name,
                       test.df = test.df)
  stopCluster(cl)
  return(results)

}

#elastic net
get_en_results <- function(sample.list, model.vars,outcome.name,test.df){
  #Fits in parallel elasticnet varaible selection to a list of dataframes
  #param: sample.list: list of dataframes containing sample data
  #param: model.vars: the variable in dataframe
  #param: outcome.name: string name of outcome
  spec_for_sens <- function(model,outcome.name,test.df){
    newx = model.matrix(formula(paste(outcome.name, "~", ".")), data = test.df)
    #compute roc
    en_prob <- predict(model, type="prob",
                       newdata  = newx , s = model$finalModel$lambdaOpt)

    r1 <- pROC::roc(test.df$final_outcome_death, en_prob$X1)  #computes ROC curve from model applied to test set
    #spec
    ss1 <- data.frame(r1$sensitivities, r1$specificities); names(ss1) <- c('Sens', "Spec")
    #sens
    spec1 <- dplyr::filter(ss1, abs(Sens - 0.8) == min(abs(Sens - 0.8)))$Spec
    #if more than 1 get mean
    if(length(spec1) > 1){
      spec1 <- mean(spec1,na.rm=T)
    }
    return(spec1)

  }
  #helper to speed up
  get_en_vars <- function(sample.df, model.vars,outcome.name,test.df){

    fit.en <- function(df,model.vars,outcome.name,test.df){


      #get n and variable names
      n <- nrow(df); df <- na.omit(df); y.idx <- grep(outcome.name,colnames(df))
      #get x and y dataframes
      df.x <- model.matrix(formula(paste(outcome.name, "~", ".")), data = df); df.y <- as.factor(df[,y.idx])
      #fit lasso
      tr.ctrl <- trainControl(method="cv", number=5,
                              classProbs=TRUE, summaryFunction=twoClassSummary,
                              allowParallel = T)

      fit.glm <- train(x = df.x,y=df.y,
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
      spec.sens <- spec_for_sens(fit.glm,outcome.name,test.df)

      return(c(in.model,"EN",n, roc.val, spec.sens))

    }

    #fit row and get vars chosen
    out.row <- as.data.frame(t(fit.en(sample.df, model.vars,outcome.name,test.df)))
    model.vars <- c(model.vars,  "method", "N", "ROC","spec_for_sens")
    names(out.row) <- model.vars
    return(out.row)

  }

  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, library(caret))
  results <- parLapply(cl, sample.list, get_en_vars,
                       model.vars = model.vars, outcome.name=outcome.name,
                       test.df=test.df)
  stopCluster(cl)

  #Fitting takes ~ 480seconds

  return(results)

}

#adapative lasso
get_adaLasso_results <- function(sample.list, model.vars,outcome.name,test.df){
  #Fits in parallel adataptive lasso to a list of dataframes
  #param: sample.list: list of dataframes containing sample data
  #param: model.vars: the variable in dataframe
  #param: outcome.name: string name of outcome
  spec_for_sens <- function(model,outcome.name,test.df){
    newx = model.matrix(formula(paste(outcome.name, "~", ".")), data = test.df)
    #compute roc
    adaLasso_prob <- predict(model, type="prob",
                             newdata  = newx , s = model$finalModel$lambdaOpt)

    r1 <- pROC::roc(test.df$final_outcome_death, adaLasso_prob$X1)  #computes ROC curve from model applied to test set
    #spec
    ss1 <- data.frame(r1$sensitivities, r1$specificities); names(ss1) <- c('Sens', "Spec")
    #sens
    spec1 <- dplyr::filter(ss1, abs(Sens - 0.8) == min(abs(Sens - 0.8)))$Spec
    #if more than 1 get mean
    if(length(spec1) > 1){
      spec1 <- mean(spec1,na.rm=T)
    }
    return(spec1)

  }

  #helper to speed up
  get_adaLasso_vars <- function(sample.df, model.vars,outcome.name,test.df){

    fit.adaLasso <- function(df,model.vars,outcome.name,test.df){


      #get n and variable names
      n <- nrow(df); df <- na.omit(df); y.idx <- grep(outcome.name,colnames(df))
      #get x and y dataframes
      df.x <- model.matrix(formula(paste(outcome.name, "~", ".")), data = df); df.y <- as.factor(df[,y.idx])
      #train control
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

      final.model <- coef(fit.adaLasso$finalModel,fit.adaLasso$bestTune$lambda)
      top.vars <- rownames(final.model)[final.model[,1]!= 0]
      #score variable with 1 if it is in the model and 0 o.w
      roc.val <- max(fit.adaLasso$results$ROC)
      in.model <- as.numeric(model.vars %in% top.vars)
      spec.sens <- spec_for_sens(fit.glm,outcome.name,test.df)
      return(c(in.model,"adaLasso",n, roc.val,spec.sens))

    }

    #fit row and get vars chosen
    out.row <- as.data.frame(t(fit.adaLasso(sample.df, model.vars,outcome.name,test.df)))
    model.vars <- c(model.vars,  "method", "N", "ROC", "spec_for_sens")
    names(out.row) <- model.vars
    return(out.row)

  }

  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, library(caret))
  results <- parLapply(cl, sample.list, get_adaLasso_vars,
                       model.vars = model.vars, outcome.name=outcome.name,
                       test.df=test.df)
  stopCluster(cl)

  #Fitting takes ~ 480seconds

  return(results)

}

get_rfe_results <- function(sample.list, model.vars,outcome.name,test.df){
  #Fits in parallel adataptive rfe to a list of dataframes
  #param: sample.list: list of dataframes containing sample data
  #param: model.vars: the variable in dataframe
  #param: outcome.name: string name of outcome
  spec_for_sens <- function(model,outcome.name,test.df){
    newx = model.matrix(formula(paste(outcome.name, "~", ".")), data = test.df)
    #compute roc
    rfe_prob <- predict(model, type="prob",
                        newdata  = newx )

    r1 <- pROC::roc(test.df$final_outcome_death, rfe_prob$X1)  #computes ROC curve from model applied to test set
    #spec
    ss1 <- data.frame(r1$sensitivities, r1$specificities); names(ss1) <- c('Sens', "Spec")
    #sens
    spec1 <- dplyr::filter(ss1, abs(Sens - 0.8) == min(abs(Sens - 0.8)))$Spec
    #if more than 1 get mean
    if(length(spec1) > 1){
      spec1 <- mean(spec1,na.rm=T)
    }
    return(spec1)

  }

  #helper to speed up
  get_rfe_vars <- function(sample.df,  model.vars,outcome.name,test.df){
    #Get the top variables chosen by RFE for different samples sizes
    #helper functions
    fit.rfe <- function(df,model.vars,outcome.name,test.df){
      #function that will fit a rfe and return top 10 most important varibales

      #get n and variable names
      n <- nrow(df); df <- na.omit(df); y.idx <- grep(outcome.name,colnames(df))
      #get x and y dataframes
      df.x <- model.matrix(formula(paste(outcome.name, "~", ".")), data = df); df.y <- as.factor(df[,y.idx])
      #fit rf
      lrFuncs$summary <- twoClassSummary;
      ctrl <- rfeControl(functions=lrFuncs, method = 'cv', repeats = 5,verbose = F, allowParallel = T,saveDetails = T)
      subsets <- seq(1,25,3)
      rfe.model <- rfe(x=df.x,y=df.y,sizes = subsets, metric = "ROC", rfeControl = ctrl)
      #get top variables
      in.model <- as.numeric(model.vars %in% rfe.model$optVariables)
      #get roc value
      idx <- rfe.model$results$Variables == rfe.model$optsize
      roc.val <- rfe.model$results[idx,2]
      spec.sens <- spec_for_sens(rfe.model,outcome.name,test.df)

      return(c(in.model,"RFE",n, roc.val,spec.sens))

    }


    #fit row and get vars chosen
    out.row <- as.data.frame(t(fit.rfe(sample.df, model.vars,outcome.name,test.df)))
    model.vars <- c(model.vars, "method", "N", "ROC", "spec_for_sens")
    names(out.row) <- model.vars
    return(out.row)

  }

  #fit in parallel and get the results
  #make cluster
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, {library(caret)})
  results <- parLapply(cl,sample.list, get_rfe_vars,
                       model.vars = model.vars, outcome.name = outcome.name,
                       test.df = test.df)
  stopCluster(cl)
  #Fitting takes ~ 480seconds

  return(results)
}
