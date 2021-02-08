#load packages 
library(tidyverse)
library(ROCR)
library(grid)
library(broom)
library(caret)
library(glmnet)
library(ROSE)
library(RColorBrewer)
library(rpart)
library(coefplot)
#library(mlr) #uncomment to attach mlr. need to dettach to use ROCR
#Smart Discharge R EDA 
df <- read_csv('data/SD_dat_for_Ash1142019.csv') 
df$final_outcome_death <- as.factor(df$final_outcome_death)
df$hivMaternal <- as.factor(df$hivMaternal) 
df <- df %>% rename(outcome = final_outcome_death)
df <- df %>% drop_na(outcome)
levels(df$outcome) <- c(0,  1)

df <- df %>% select(-c(dpi,Dataset))
#######################################################################################################################
#table of outcomes
table(df$outcome) #big imbalance between classes. Approx ~5% of outcomes are death 
#Jeff suggestions/comments
#1.	Create a single tree model which could be used to reasonably predict death from these set of variables.
#The lower number of levels the better. Use AUROC as your model choosing criteria, 
#because if you use ‘accuracy’ is will be based off a 50% cutoff for risk, and no patients risk will be that high
#
#2.	Determine for each variable in the dataset if and what a cutoff value would be that could reasonably stratify risk.
#Model building:
#Before constructing a tree based model try to get the best logistic regression model possible
#Then compare the preformance of tree based methods to this simple model. 
#Methods:
#1. Construct simple logistic regression
#2. Logistic regression with some form of variable selection 
#3. Logistc regreession with some penalty function 
#4. Logistic regression with over/undersampling methods
#Due to the dataset being unbalanced, a simple p>0.5 to classify will not work very well.
#Need to determine and optimal cut off point for any logistic regression model that is fit
#precsion -> class is correctly recognised (tp / tp + fp7/)
#recall -> postive actually postive (tp/ sum(p))
#High recall, low precsion -> lots of postives correctly recongnized (low fn) but lots of fp
#low recall, high precsion -> high fn, but low fp

#######################################################################################################################
#1.Logistic regression model
cv_logistic_cutOff <- function(cutOff, model, k, df){
  #function that fits a logistic regression model for a given cuttOff
  #using cross validation to determine the precsion and recall over 
  #k-folds of cross validation 
  #@param cutOff : cuttOff value
  #@param model: model 
  #@param k: number of folds
  #@param df: a dataframe with the outcome column named outcome


  cm_list <- list(k)
  precsion_vec <- vector(length=k)
  recall_vec <- vector(length=k)
  fp_ratio_vec <- vector(length=k)
  auc.score <- vector(length = k)
  flds <- createFolds(df$outcome, k = k, list = T, returnTrain = F)
  for(i in 1:k){
    
    #create a train and test set
    test.idx <- flds[[i]]
    df.train <- df[-test.idx,]
    df.test <- df[test.idx,]
    #fit model
    fit.logit <- glm(model, data = df.train, family = binomial(link='logit'))
    #predictions
    y.test <- df.test$outcome
    x.test <- df.test %>% select(-outcome)
    pred.prob <- predict(fit.logit,x.test, type='response')
    predictions <- rep(0,length(pred.prob))
    idx_1 <- pred.prob > cutOff
    predictions[idx_1] = 1  
    predictions <- as.factor(predictions)
    #get auc 
    pred.auc <- prediction(as.numeric(predictions), as.numeric(y.test))
    auc.val <- performance(pred.auc, 'auc')
    auc.score[i] <- as.numeric(slot(auc.val, 'y.values'))
    #get confusion matrix
    cm <- confusionMatrix(predictions, y.test, positive = '1', mode='everything')
    cm_list[[i]] <- cm
    precsion_vec[i] <- cm$byClass[5]
    recall_vec[i] <- cm$byClass[6]
    n <- cm$table[1,1]
    fp <- cm$table[2,1]
    fp_ratio_vec[i] <- fp / (fp + n)
    
  }
  
  mean_precsion <- mean(precsion_vec)
  mean_recall <- mean(recall_vec)
  mean_fp_ratio <- mean(fp_ratio_vec)
  mean_auc <- mean(auc.score)
  
  return(c(mean_precsion,mean_recall,mean_fp_ratio, mean_auc))
  
}

plot_model_cuttOff <- function(model,k,title, grid,df){
  #function to cross validate a grid of cuttoff values for a given model
  grid_cutoff <- grid
  output_cv <- sapply(grid_cutoff, cv_logistic_cutOff, model = model, k = k, df=df)
  precion_out <- output_cv[1,]
  recall_out <- output_cv[2,]
  fp_ratio <- output_cv[3,]
  auc_score <- output_cv[4,]
  plot_df <- tibble(precsion = precion_out, recall=recall_out,grid_vals = grid_cutoff, fp_ratio = fp_ratio, auc = auc_score)
  plot_df %>% 
    ggplot(aes(grid_cutoff))+
    geom_line(aes(y = precsion, colour = 'precsion'))+
    geom_line(aes(y = recall, colour='recall'))+
    geom_line(aes(y = fp_ratio, colour = 'fp_ratio'))+
    geom_line(aes(y = auc, colour='auc'))+
    theme_bw()+
    ggtitle(title)
  
}
#######################################################################################################################
#rename df for functions
k <- 10
#Good candiates from previous work: Distancenew, W_AZ, BMI_AZ, W_LZ
p_distancenew <- plot_model_cuttOff(as.formula("outcome ~ Distancenew"),10, "Distancenew only ", grid = seq(0,0.2,0.005), df)
p_w_az <- plot_model_cuttOff(as.formula("outcome ~ W_AZ"),10, "W_AZ only ", grid = seq(0,0.15,0.005), df)
p_bmi_az <- plot_model_cuttOff(as.formula("outcome ~ BMI_AZ"),10, "BMI_AZ only ", grid = seq(0,0.15,0.005), df)
p_w_lz <- plot_model_cuttOff(as.formula("outcome ~ W_LZ"),10, "W_LZ only ", grid = seq(0,0.15,0.005), df)
p_muac <- plot_model_cuttOff(as.formula("outcome ~ muac"),10, "muac only ", grid = seq(0,0.15,0.005), df)
p_hg <-  plot_model_cuttOff(as.formula("outcome ~ hg"),10, "hg only ", grid = seq(0,0.15,0.005), df)
p_distancenew;p_w_az;p_bmi_az;p_w_lz;p_muac;p_hg
#distance new seems to be fairly decent 

p_full <-  plot_model_cuttOff(as.formula("outcome ~ ."),10, "full mode", grid = seq(0,0.5,0.005), df)
p_full

#using only distance, W_AZ, muac and hg
p_lasso <-  plot_model_cuttOff(as.formula("outcome ~ ."),10, "LASSO mode", grid = seq(0,0.5,0.005), 
                              df %>% select(outcome, muac, W_AZ, hg))#df.small - only LASSO vars
p_lasso
#######################################################################################################################
#fitting a penalized logistic regresion model 
#running 10 fold cv
logistic_en.cutOff <- function(cutOff, k, df,alpha.val){
  #function that runs k-cold cv and fits an elastic net logisic regression.
  #then determines predicted class for some given cutOff value
  #need to supply cutOff value, k, a df and the alpha value for elastic net
  
  #create folds
  set.seed(10)
  flds <- createFolds(df$outcome, k = k, list = T, returnTrain = F)
  precsion_vec <- vector(length=k)
  recall_vec <- vector(length=k)
  fp_ratio_vec <- vector(length=k)
  for(i in 1:k){
    
    #create a train and test set
    test.idx <- flds[[i]]
    df.train <- df[-test.idx,]
    x.train <- df.train %>% select(-outcome) %>% sapply(., as.numeric)
    y.train <- df.train %>% select(outcome) 
    y.train <- y.train$outcome
    df.test <- df[test.idx,]
    x.test <- df.test %>% select(-outcome) %>% sapply(., as.numeric)
    y.test <- df.test %>% select(outcome) 
    y.test <- y.test$outcome
    
    #fit elastic net
    lambdas <- exp(seq(-3,10,length=50))
    fit.en <- cv.glmnet(x.train,y.train,family='binomial', lambda = lambdas, nfolds = 10,alpha=alpha.val)
    
    #get predicted probabilties
    pred.prob <- predict(fit.en,x.test,s="lambda.min" ,type='response')
    
    #get predictions given cutOff
    pred <- rep(0,length(pred.prob))
    idx_1 <- pred.prob > cutOff
    pred[idx_1] = 1
    
    #make confusion matrix
    cm <- confusionMatrix(factor(pred),factor(y.test),positive = '1', mode='everything')
    precsion_vec[i] <- cm$byClass[5]
    recall_vec[i] <- cm$byClass[6]
    tn <- cm$table[1,1]
    fp <- cm$table[2,1]
    fp_ratio_vec[i] <- fp / (fp + tn)
    
  }
  
  mean_precsion <- mean(precsion_vec)
  mean_recall <- mean(recall_vec)
  mean_fp_ratio <- mean(fp_ratio_vec)
  
  return(c(mean_precsion,mean_recall,mean_fp_ratio))

}

plot_en_model_cuttOff <- function(output_cv,grid_cutoff,title){
  #function to cross validate a grid of cuttoff values for a given model
  precion_out <- output_cv$precsion
  recall_out <- output_cv$recall
  fp_ratio <- output_cv$fp_ratio
  plot_df <- tibble(precsion = precion_out, recall=recall_out,grid_vals = grid_cutoff, fp_ratio = fp_ratio)
  plot_df %>% 
    ggplot(aes(grid_cutoff))+
    geom_line(aes(y = precsion, colour = 'precsion'))+
    geom_line(aes(y = recall, colour='recall'))+
    geom_line(aes(y = fp_ratio, colour = 'fp_ratio'))+
    theme_bw()+
    ggtitle(title)
  
}

#######################################################################################################################
k <- 10
grid_cutOff <- seq(0,0.15,0.005)


output_.05 <- sapply(grid_cutOff, logistic_en.cutOff, k = k, df = df, alpha.val = 0.05)
output_.05 <- t(output_.05)
output_.05 <- tibble(grid = grid_cutOff, precsion = output_.05[,1], fp_ratio= output_.05[,3], recall = output_.05[,2] )
plot_en_model_cuttOff(output_cv = output_.05,grid_cutoff = grid_cutOff,title = 'EN: alpha = 0.05')

output_.25 <- sapply(grid_cutOff, logistic_en.cutOff, k = k, df = df, alpha.val = 0.25)
output_.25 <- t(output_.25)
output_.25 <- tibble(grid = grid_cutOff, precsion = output_.25[,1], fp_ratio= output_.25[,3], recall = output_.25[,2] )
plot_en_model_cuttOff(output_cv = output_.25,grid_cutoff = grid_cutOff,title = 'EN: alpha = 0.25')


output_.95 <- sapply(grid_cutOff, logistic_en.cutOff, k = k, df = df, alpha.val = 0.95)
output_.95 <- t(output_.95)
output_.95 <- tibble(grid = grid_cutOff, precsion = output_.95[,1], fp_ratio= output_.95[,3], recall = output_.95[,2] )
plot_en_model_cuttOff(output_cv = output_.95,grid_cutoff = grid_cutOff,title = 'EN: alpha = 0.95')


output_.75 <- sapply(grid_cutOff, logistic_en.cutOff, k = k, df = df, alpha.val = 0.75)
output_.75 <- t(output_.75)
output_.75 <- tibble(grid = grid_cutOff, precsion = output_.75[,1], fp_ratio= output_.75[,3], recall = output_.75[,2] )
plot_en_model_cuttOff(output_cv = output_.75,grid_cutoff = grid_cutOff,title = 'EN: alpha = 0.75')


output_.5 <- sapply(grid_cutOff, logistic_en.cutOff, k = k, df = df, alpha.val = 0.50)
output_.5 <- t(output_.5)
output_.5 <- tibble(grid = grid_cutOff, precsion = output_.5[,1], fp_ratio= output_.5[,3], recall = output_.5[,2] )
plot_en_model_cuttOff(output_cv = output_.5,grid_cutoff = grid_cutOff,title = 'EN: alpha = 0.50')

output_ridge <- sapply(grid_cutOff, logistic_en.cutOff, k = k, df = df, alpha.val = 0)
output_ridge <- t(output_ridge)
output_ridge <- tibble(grid = grid_cutOff, precsion = output_ridge[,1], fp_ratio= output_ridge[,3], recall = output_ridge[,2] )
plot_en_model_cuttOff(output_cv = output_ridge,grid_cutoff = grid_cutOff,title = 'Ridge regression')


#######################################################################################################################
#Thus far: 
#Logistic regression has fairly poor preformance. Elastic Net approach wasn't much better. Ridge seems to preform better than LASSO.

#Next: 
#- Random oversampling: Take a random sample of minority group with replacement to increase miniority group size


#Oversampling 
df.overSamp <- ovun.sample(outcome ~ ., data = df,method='over', p = 0.4)$data
#fit logistic regression on this data
k <- 10
#Good candiates from previous work: Distancenew, W_AZ, BMI_AZ, W_LZ
p_distancenew <- plot_model_cuttOff(as.formula("outcome ~ Distancenew"),10, "Distancenew only ", grid = seq(0,0.8,0.05), df.overSamp)
p_w_az <- plot_model_cuttOff(as.formula("outcome ~ W_AZ"),10, "W_AZ only ", grid = seq(0,0.8,0.05), df.overSamp)
p_bmi_az <- plot_model_cuttOff(as.formula("outcome ~ BMI_AZ"),10, "BMI_AZ only ", grid = seq(0,0.8,0.05), df.overSamp)
p_w_lz <- plot_model_cuttOff(as.formula("outcome ~ W_LZ"),10, "W_LZ only ", grid = seq(0,0.8,0.05), df.overSamp)
p_muac <- plot_model_cuttOff(as.formula("outcome ~ muac"),10, "muac only ", grid = seq(0,0.8,0.05), df.overSamp)
p_distancenew;p_w_az;p_bmi_az;p_w_lz;p_muac
#distance new seems to be fairly decent 
p_full <-  plot_model_cuttOff(as.formula("outcome ~ ."),10, "full mode", grid = seq(0,0.8,0.05), df.overSamp)
p_full


#Oversampling 
df.overSamp <- ovun.sample(outcome ~ ., data = df,method='over')$data
#fit logistic regression on this data
k <- 10
#Good candiates from previous work: Distancenew, W_AZ, BMI_AZ, W_LZ
p_distancenew <- plot_model_cuttOff(as.formula("outcome ~ Distancenew"),10, "Distancenew only ", grid = seq(0,0.8,0.05), df.overSamp)
p_w_az <- plot_model_cuttOff(as.formula("outcome ~ W_AZ"),10, "W_AZ only ", grid = seq(0,0.8,0.05), df.overSamp)
p_bmi_az <- plot_model_cuttOff(as.formula("outcome ~ BMI_AZ"),10, "BMI_AZ only ", grid = seq(0,0.8,0.05), df.overSamp)
p_w_lz <- plot_model_cuttOff(as.formula("outcome ~ W_LZ"),10, "W_LZ only ", grid = seq(0,0.8,0.05), df.overSamp)
p_muac <- plot_model_cuttOff(as.formula("outcome ~ muac"),10, "muac only ", grid = seq(0,0.8,0.05), df.overSamp)
p_distancenew;p_w_az;p_bmi_az;p_w_lz;p_muac
#distance new seems to be fairly decent 
p_full <-  plot_model_cuttOff(as.formula("outcome ~ ."),10, "full mode", grid = seq(0,0.8,0.05), df.overSamp)
p_full

#ROSE method
df.overSamp <- ROSE(outcome ~ ., data = df,seed=1)$data
#fit logistic regression on this data
k <- 10
#Good candiates from previous work: Distancenew, W_AZ, BMI_AZ, W_LZ
p_distancenew <- plot_model_cuttOff(as.formula("outcome ~ Distancenew"),10, "Distancenew only ", grid = seq(0,0.8,0.05), df.overSamp)
p_w_az <- plot_model_cuttOff(as.formula("outcome ~ W_AZ"),10, "W_AZ only ", grid = seq(0,0.8,0.05), df.overSamp)
p_bmi_az <- plot_model_cuttOff(as.formula("outcome ~ BMI_AZ"),10, "BMI_AZ only ", grid = seq(0,0.8,0.05), df.overSamp)
p_w_lz <- plot_model_cuttOff(as.formula("outcome ~ W_LZ"),10, "W_LZ only ", grid = seq(0,0.8,0.05), df.overSamp)
p_muac <- plot_model_cuttOff(as.formula("outcome ~ muac"),10, "muac only ", grid = seq(0,0.8,0.05), df.overSamp)
p_distancenew;p_w_az;p_bmi_az;p_w_lz;p_muac
#distance new seems to be fairly decent 
p_full <-  plot_model_cuttOff(as.formula("outcome ~ ."),10, "full mode", grid = seq(0,0.8,0.05), df.overSamp)
p_full


#######################################################################################################################
#oversampling and then elastic net
df.overSamp <- ovun.sample(outcome ~ ., data = df,method='over')$data

k <- 10
grid_cutOff <- seq(0,0.8,0.05)

output_.25 <- sapply(grid_cutOff, logistic_en.cutOff, k = k, df = df.overSamp, alpha.val = 0.25)
output_.25 <- t(output_.25)
output_.25 <- tibble(grid = grid_cutOff, precsion = output_.25[,1], fp_ratio= output_.25[,3], recall = output_.25[,2] )
plot_en_model_cuttOff(output_cv = output_.25,grid_cutoff = grid_cutOff,title = 'EN: alpha = 0.25')


output_.75 <- sapply(grid_cutOff, logistic_en.cutOff, k = k, df = df.overSamp, alpha.val = 0.75)
output_.75 <- t(output_.75)
output_.75 <- tibble(grid = grid_cutOff, precsion = output_.75[,1], fp_ratio= output_.75[,3], recall = output_.75[,2] )
plot_en_model_cuttOff(output_cv = output_.75,grid_cutoff = grid_cutOff,title = 'EN: alpha = 0.75')


output_.5 <- sapply(grid_cutOff, logistic_en.cutOff, k = k, df = df.overSamp, alpha.val = 0.50)
output_.5 <- t(output_.5)
output_.5 <- tibble(grid = grid_cutOff, precsion = output_.5[,1], fp_ratio= output_.5[,3], recall = output_.5[,2] )
plot_en_model_cuttOff(output_cv = output_.5,grid_cutoff = grid_cutOff,title = 'EN: alpha = 0.50')

output_ridge <- sapply(grid_cutOff, logistic_en.cutOff, k = k, df = df.overSamp, alpha.val = 0)
output_ridge <- t(output_ridge)
output_ridge <- tibble(grid = grid_cutOff, precsion = output_ridge[,1], fp_ratio= output_ridge[,3], recall = output_ridge[,2] )
plot_en_model_cuttOff(output_cv = output_ridge,grid_cutoff = grid_cutOff,title = 'Ridge regression')

#######################################################################################################################
#[Update]
#Ridge regression seems to be the best. 
#Need to cross validate and compare with simple logistic regression 

get_recall <- function(cm){
  #takes cm returns recall
  #         |predicted|
  #         |0 | 1 |
  #actual  0|__|__ |
  #        1|fn|tp |
  
  
  tp <- cm[2,2]
  fp <- cm[2,1]
  
  recall <- tp / (tp + fp)
  
  return(recall)
}

get_frat <- function(cm){
  #takes cm returns fp_rate
  #         |predicted|
  #         |0 | 1 |
  #actual  0|tn|fp |
  #        1|  |   |
  
  fp <- cm[1,2]
  tn <- cm[1,1]
  
  frat <- fp / (fp + tn)
  
  return(frat)
  
}


N <- 50
#containers for mean CV recall and frat
pred.ridge.c1.recall <- pred.ridge.c2.recall <-pred.log.c1.recall <- pred.log.c2.recall <- rep(0,N)
pred.ridge.c1.frat <- pred.ridge.c2.frat <-pred.log.c1.frat <- pred.log.c2.frat <- rep(0,N)
for(j in 1:N){
  #Do N-runs of 10-fold CV
  
  #declare k-folds
  k <- 10
  #make folds
  flds <- createFolds(df$outcome, k = 10, list = T, returnTrain = F)
  #container objects for predictions
  pred.ridge.c1 <- list(10)
  pred.ridge.c2 <- list(10)
  pred.log.c1 <- list(10)
  pred.log.c2 <- list(10)
  
  for(i in 1:k){
    
    #create a train and test set 
    test.idx <- flds[[i]]
    df.train <- df[-test.idx,]
    df.test <- df[test.idx,]
    
    #fit logistic regression model
    x.test.logistic <- df.test %>% select(-outcome)
    fit.logistic <- glm(outcome ~ .,data=df.train, family=binomial(link = 'logit'))
    pred.prob.logistic <- predict(fit.logistic, x.test.logistic,type='response')
    
    #fit ridge regression
    x.train.ridge <-  df.train %>% select(-outcome) %>% sapply(., as.numeric)
    x.test.ridge <- df.test %>% select(-outcome) %>% sapply(., as.numeric)
    fit.ridge <- cv.glmnet(x.train.ridge,df.train$outcome, family='binomial',
                           lambda = exp(seq(-3,10,length=50)), nfolds = 10, alpha=0)
    pred.prob.ridge <- predict(fit.ridge,x.test.ridge, s = 'lambda.min', type='response')
    
    
    #get predictions
    pred.vec1 <- pred.vec2 <- pred.vec3 <- pred.vec4 <-  rep(0,nrow(df.test))
    
    #use optCutoffs found previously
    pred.vec1[pred.prob.ridge > 0.04] <- 1
    pred.vec2[pred.prob.ridge > 0.045] <- 1
    pred.vec3[pred.prob.logistic > 0.025] <- 1
    pred.vec4[pred.prob.logistic > 0.035] <- 1
    
    #save predictions
    pred.ridge.c1[[i]] <-  table(df.test$outcome, pred.vec1)
    pred.ridge.c2[[i]] <-  table(df.test$outcome, pred.vec2)
    pred.log.c1[[i]] <-  table(df.test$outcome, pred.vec3)
    pred.log.c2[[i]] <-  table(df.test$outcome, pred.vec4)
    
      
  }
  
  
  #get recall and frat for all k-CV runs
  pred.ridge.c1.recall[j] <- mean(sapply(pred.ridge.c1, get_recall))
  pred.ridge.c2.recall[j] <- mean(sapply(pred.ridge.c2, get_recall))
  pred.log.c1.recall[j] <- mean(sapply(pred.log.c1, get_recall))
  pred.log.c2.recall[j] <- mean(sapply(pred.log.c2, get_recall))
  
  pred.ridge.c1.frat[j] <- mean(sapply(pred.ridge.c1, get_frat))
  pred.ridge.c2.frat[j] <- mean(sapply(pred.ridge.c2, get_frat))
  pred.log.c1.frat[j] <- mean(sapply(pred.log.c1, get_frat))
  pred.log.c2.frat[j] <- mean(sapply(pred.log.c2, get_frat))
  
  
}


#Plot histograms of recall and frat 
recall_df <- tibble(logc1 = pred.log.c1.recall, logc2 = pred.log.c2.recall,
                    ridgec1 = pred.ridge.c1.recall, ridgec2 = pred.ridge.c2.recall)

frat_df <- tibble(logc1 = pred.log.c1.frat, logc2 = pred.log.c2.frat,
                  ridgec1 = pred.ridge.c1.frat, ridgec2 = pred.ridge.c2.frat)

recall_df %>% 
  gather(.) %>% 
  ggplot(., aes(x = key, y = value, fill = factor(key) ))+
  geom_boxplot(width=0.4)+
  ylab('Recall')+
  ggtitle('Recall for different regression models')+
  theme_bw()+
  scale_fill_brewer(palette = 'Dark2')


frat_df %>% 
  gather(.) %>% 
  ggplot(., aes(x = key, y = value, fill = factor(key) ))+
  geom_boxplot(width=0.4)+
  ylab('fp_ratio')+
  ggtitle('fp_ratio for different regression models')+
  theme_bw()+
  scale_fill_brewer(palette = 'Dark2')

#######################################################################################################################
#Now have a standard to compare to.
#Note: Confusion Matrix vs AUC evaluation of fit
#A confusion matrix evaluates one particular classifier with a fixed threshold, while the AUC evaluates that classifier over all possible thresholds.
#Going to use the above metric to compare model preformance. 
#######################################################################################################################
#Decsion Trees: 
#Logistic regression was fairly poor at classifying. Maybe a decsion tree can do better
#Resources:
#https://rstudio-pubs-static.s3.amazonaws.com/442284_82321e66af4e49d58adcd897e00bf495.html
#https://cran.r-project.org/web/packages/rpart/vignettes/longintro.pdf
#http://rstudio-pubs-static.s3.amazonaws.com/3588_81e2ebd4de1b41bc9ac2f29f5f7dab2e.html
#https://www.gormanalysis.com/blog/decision-trees-in-r-using-rpart/
#CART Trees:
#1. Build basic CART Tree (using different node purity measures)
#2. Build CART tree by providing a loss matrix that changes the
#   missclassification penalty, such that a fn is penalized heavily.
#3. Try ADABoosting/other resampling/weighting methods
#4. PRIM/PolyMARS methods instead of classic CART
#5. Try HME if there is time

#make default tree
fit.info <- rpart(outcome ~ . , data = df,method = 'class', parms = list(split = 'information'))
fit.gini <- rpart(outcome ~ ., data = df, method = 'class',parms = list(split = 'gini'))
#rpart was unable to fit a tree. Need to change control parameters
control.deafult <- rpart.control(minsplit = 5, minbucket = 5, cp = 0)
fit2.info <- rpart(outcome ~ . , data = df,method = 'class', parms = list(split = 'information'), control = "")
fit2.gini <- rpart(outcome ~ . , data = df,method = 'class', parms = list(split = 'gini'), control = control.deafult)
#generate confusion matrix on dataset to assses gof
info.pred <- predict(fit2.info,type='class')
gini.pred <- predict(fit2.gini, type='class')
table(info.pred, df$outcome); table(gini.pred,df$outcome)
#gini appears to have a higher rate of FN but fewer FP. 
#info gain had a lower rate of FN, which is what we care out
#Note this is just prelim, without any cross validation work.


#make tree with different loss matricies
loss.1 <- matrix(c(0,20,1,0), nrow=2)
fit1.loss <- rpart(outcome ~ ., data = df, method = 'class', parms = list(split = 'information', loss =loss.1 ))
fit2.loss <- rpart(outcome ~ ., data = df, method = 'class', parms = list(split = 'gini', loss =loss.1 ))
#generate confusion matrix on dataset to assses gof
loss1.pred <- predict(fit1.loss,type='class')
loss2.pred <- predict(fit2.loss, type='class')
table(loss1.pred, df$outcome); table(loss2.pred, df$outcome)
#get auc
prediction.fit1 <- prediction(as.numeric(loss1.pred), as.numeric(df$outcome))
auc.fit1 = performance(prediction.fit1,'auc')
#add a loss matrix really improves preformance. fn ~ 25%, down from 37%
#changing loss matrix to: loss.1 <- matrix(c(0,20,1,0), nrow=2), gives 0% fn
#NB: (might be massively overfitting) --> need to evaluate using CV
#######################################################################################################################
#DT not great. Could try with fewer variables
#Var selection;
#fit elastic net
lambdas <- exp(seq(-20,10,length=500))
x.train <- model.matrix( ~ ., dplyr::select(df, -outcome))[, -1]
y.train <- df %>% select(outcome)
fit.en <- cv.glmnet(x.train,y.train$outcome,family='binomial', lambda = lambdas,nfolds = 10,alpha=1)
plot(fit.en, )
#get coefs of sparse model
lambda.min <- fit.en$lambda.min
lasso_model <- extract.coef(fit.en, s = lambda.min)
lasso_model <- as.character(lasso_model$Coefficient)[-1]
lasso_model[10] <- "hivMaternal"

#function to CV a decsion tree
cv_classTree <- function(df,k, parms, control, formula){
  
  flds <- createFolds(df$outcome, k = k, list = T, returnTrain = F)
  precsion_vec <- recall_vec <- fpRat_vec <- auc_vec <-  rep(0,k)
  for(i in 1:k){
    #create a train and test set
    test.idx <- flds[[i]]
    df.train <- df[-test.idx,]
    df.test <- df[test.idx,]
    #get test x and y
    y.test <- df.test %>% select(outcome)
    x.test <- df.test %>% select(-outcome)
    #fit tree
    tree.fit <- rpart(formula = formula, data = df.train, method = 'class', 
                      parms = parms, control = control)
    cp.min <- tree.fit$cptable[,1][which.min(tree.fit$cptable[,4])]
    tree.pruned <- prune.rpart(tree.fit, cp.min)
    #make predictions on test set
    tree.pred <- predict(tree.pruned,newdata = x.test , type = 'class')
    #get AUC for decsion tree
    rocr.prediction <- prediction(as.numeric(tree.pred), as.numeric(y.test$outcome))
    auc.pred <- performance(rocr.prediction, 'auc')
    auc_vec[i] <- as.numeric(slot(auc.pred, 'y.values'))
    #generate confusion matrix and asses prediction accuracy
    cm <- table(tree.pred, y.test$outcome)
    tn <-cm[1,1]; fn <- cm[1,2]; fp <- cm[2,1]; tp <- cm[2,2]
    precsion <- tp/(tp + fp); recall <- tp/(tp + fn); 
    fp_rat <- fp/(fp + tn)
    precsion_vec[i] <- precsion; recall_vec[i] <- recall
    fpRat_vec[i] <- fp_rat
    
  }
  
  return_df <- tibble(precsion = precsion_vec, recall = recall_vec, fp_rat = fpRat_vec, auc = auc_vec)
  return(return_df)
  
}
df.small <- df %>%  select(outcome, muac, W_AZ, hg)
set.seed(188)



K <- 10
cv.precsion <- cv.recall <- cv.fp_rat <- cv.auc <- rep(0,K)
for(j in 1:K){
  
  #parameters for rpart
  loss.1 <- matrix(c(0,10^10,1,0), nrow=2)
  k <- 10
  parms <- list(split = 'information', loss =loss.1 )
  formula.fit <- as.formula("outcome ~ .")
  control <- rpart.control(minsplit = 25,minbucket = 10,cp = 0.02,maxdepth = 20)
  
  
  #fit tree
  cv_out <- cv_classTree(df = df.small, k = k, parms = parms, formula = formula.fit, control = control)
  
  #get preformance metrics
  mean_cv <- cv_out %>% summarise_all(., mean)
  cv.precsion[j] <- mean_cv$precsion; cv.recall[j] <- mean_cv$recall
  cv.fp_rat[j] <- mean_cv$fp_rat; cv.auc[j] <- mean_cv$auc
}


#plot preformance results
cv_df <- tibble(cv.precsion = cv.precsion,
                cv.recall = cv.recall,
                cv.fr_rat = cv.fp_rat,
                cv.auc = cv.auc)
cv_df %>% 
  gather(.) %>% 
  ggplot(.,aes(x =key, y = value, fill= factor(key)))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_brewer(palette = 'Dark2')+
  facet_wrap( key ~ .,scales = "free")


cv_df %>% summarise_all(mean)

#######################################################################################################################
#Try using mlr package to search grid of values for optimal parameters
idx.tr <- sample(1:nrow(df.small), size = nrow(df.small)*0.85)
train  <- df.small[idx.tr,]
test <- df.small[-idx.tr,]
mort.task <- makeClassifTask(id = 'mort-train', data = train, target = 'outcome')
test.task <- makeClassifTask(id = 'mort-test', data = test, target = 'outcome')
lrn <- makeLearner('classif.rpart')  
ps <- makeParamSet(
  makeNumericParam("cp", lower = 0, upper = 0.1),
  makeNumericParam("minsplit", lower = 5, upper = 50),
  makeNumericParam("minbucket", lower = 5, upper = 20),
  makeNumericParam("maxdepth", lower = 20, upper = 30)
)

ctrl.grid <- makeTuneControlGrid()
resamp.cv <- makeResampleDesc("CV", iters= 10)
set.seed(18)
tuned.mort <- tuneParams(lrn, task = mort.task, resampling = resamp.cv,
                         control = ctrl.grid, par.set = ps,measures = list(tpr, fnr))

opt.grid <- as.data.frame(tuned.mort$opt.path)

tuned.mort2 <- tuneParams(lrn, task = mort.task, resampling = resamp.cv,
                         control = ctrl.grid, par.set = ps,measures = list(fn))

#######################################################################################################################
#Best DT was able to get: auc ~ -.65, recall ~ 0.73, fp_rat ~ 0.45






