#Smart Discharge 
#Script which will preform univariate logistic regression 
#Then fit a decsion tree 

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
library(randomForest)
########################################################################################################################################
#Functions.R
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
  specificity.vec <- vector(length = k)
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
    specificity.vec[i] <- cm$byClass[2]
    
  }
  
  mean_precsion <- mean(precsion_vec)
  mean_recall <- mean(recall_vec)
  mean_fp_ratio <- mean(fp_ratio_vec)
  mean_auc <- mean(auc.score)
  mean_spec <- mean(specificity.vec)
  
  return(c(mean_precsion,mean_recall,mean_fp_ratio, mean_auc, mean_spec))
  
}

plot_model_cuttOff <- function(model,k,title, grid,df){
  #function to cross validate a grid of cuttoff values for a given model
  grid_cutoff <- grid
  output_cv <- sapply(grid_cutoff, cv_logistic_cutOff, model = model, k = k, df=df)
  precion_out <- output_cv[1,]
  recall_out <- output_cv[2,]
  fp_ratio <- output_cv[3,]
  auc_score <- output_cv[4,]
  spec_score <- output_cv[5,]
  plot_df <- tibble(precsion = precion_out, recall=recall_out,
                    grid_vals = grid_cutoff, fp_ratio = fp_ratio, 
                    auc = auc_score, spec = spec_score)
  plot_df %>% 
    ggplot(aes(grid_cutoff))+
    geom_line(aes(y = precsion, colour = 'precsion'))+
    geom_line(aes(y = recall, colour='recall'))+
    geom_line(aes(y = fp_ratio, colour = 'fp_ratio'))+
    geom_line(aes(y = auc, colour='auc'))+
    geom_line(aes(y = spec, colour='spec'))+
    theme_bw()+
    ggtitle(title)
  
}

cv_classTree <- function(df,k, parms, control, formula){
  
  flds <- createFolds(df$outcome, k = k, list = T, returnTrain = F)
  precsion_vec <- recall_vec <- fpRat_vec <- auc_vec <- spec_vec <-  rep(0,k)
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
    #get spec
    dat.labels <- y.test$outcome == 0
    spec_vec[i] <- tn/sum(dat.labels)
    
  }
  
  return_df <- tibble(precsion = precsion_vec, recall = recall_vec, fp_rat = fpRat_vec, auc = auc_vec, spec = spec_vec)
  return(return_df)
  
}
########################################################################################################################################
#load and clean up data
df <- read_csv('data/SD_dat_for_Ash1142019.csv') 
df$final_outcome_death <- as.factor(df$final_outcome_death)
df$hivMaternal <- as.factor(df$hivMaternal) 
df <- df %>% rename(outcome = final_outcome_death)
df <- df %>% drop_na(outcome)
levels(df$outcome) <- c(0,  1)
df <- df %>% select(-c(dpi,Dataset))
df$sex <- as.factor(df$sex)
df$Abnormal_BCS <- as.factor(df$Abnormal_BCS)
df$sibling_deaths <- as.factor(df$sibling_deaths)
df$boiling <- as.factor(df$boiling)
df$hivFinal <- as.factor(df$hivFinal)
df$hivMaternal <- as.factor(df$hivMaternal)
df$bednet <- as.factor(df$bednet)
df$malaria <- as.factor(df$malaria)
#df$maternal_educationnew <- as.factor(df$maternal_educationnew)
#df$Distancenew <- as.factor(df$Distancenew)
#df$water_source_drinking <- as.factor(df$water_source_drinking)

#table of outcomes
table(df$outcome) #big imbalance between classes. Approx ~5% of outcomes are death 
#############################################################################################################################

#univariate logistic regression to determine optimal cutoffs 
k <- 10
#Good candiates from previous work: Distancenew, W_AZ, BMI_AZ, W_LZ, hg
p_distancenew <- plot_model_cuttOff(as.formula("outcome ~ Distancenew"),10, "Distancenew only ", grid = seq(0,0.05,0.0005), df)
p_distancenew
p_distancenew.opt <- plot_model_cuttOff(as.formula("outcome ~ Distancenew"),10, "Distancenew only ", grid = seq(0.03,0.04,0.0001), df)
p_distancenew.opt


p_w_az <- plot_model_cuttOff(as.formula("outcome ~ W_AZ"),10, "W_AZ only ", grid = seq(0,0.15,0.005), df)
p_w_az
p_w_az.opt <- plot_model_cuttOff(as.formula("outcome ~ W_AZ"),10, "W_AZ only ",grid = seq(0.03,0.04,0.0001), df)
p_w_az.opt

p_bmi_az <- plot_model_cuttOff(as.formula("outcome ~ BMI_AZ"),10, "BMI_AZ only ", grid = seq(0,0.15,0.005), df)
p_bmi_az

p_w_lz <- plot_model_cuttOff(as.formula("outcome ~ W_LZ"),10, "W_LZ only ", grid = seq(0,0.15,0.005), df)
p_w_lz

p_muac <- plot_model_cuttOff(as.formula("outcome ~ muac"),10, "muac only ", grid = seq(0,0.15,0.005), df)
p_muac
p_muac.opt <- plot_model_cuttOff(as.formula("outcome ~ muac"),10, "muac only ", grid = seq(0.04,0.045,0.0001), df)
p_muac.opt

p_maternaleducation <- plot_model_cuttOff(as.formula("outcome ~ maternal_educationnew"),10,
                                          "maternal_educationnew only ", grid = seq(0,0.05,0.0005), df)
p_maternaleducation


p_hg <-  plot_model_cuttOff(as.formula("outcome ~ hg"),10, "hg only ", grid = seq(0,0.15,0.005), df)
p_hg


#multivariate analysis
#fit lasso to preform model selection 
lambdas <- exp(seq(-20,10,length=500))
df.train <- model.matrix( ~ ., df)
x.train <- df.train[,-c(1,30)]
y.train <- df.train[,30]
fit.en <- cv.glmnet(x.train,y.train,family='binomial', lambda = lambdas,nfolds = 10,alpha=1)
plot(fit.en, )
#get coefs of sparse model
lambda.min <- fit.en$lambda.min
lasso_model <- extract.coef(fit.en, s = lambda.min)
lasso_model <- as.character(lasso_model$Coefficient)[-1]
lasso_model

p_full <-  plot_model_cuttOff(as.formula("outcome ~ ."),10, "Full model", grid = seq(0.01,0.025,0.0001), df)
p_full

df.lasso <- df %>% select(sex)
p_lasso <-  plot_model_cuttOff(as.formula("outcome ~ ."),10, "Model with vars selected from lasso", grid = seq(0.01,0.04,0.0005), df.lasso)
p_lasso

#using only distance, W_AZ, muac and hg
p_simple <-  plot_model_cuttOff(as.formula("outcome ~ ."),10, "Simple model", grid = seq(0.02,0.03,0.0001), 
                               df %>% select(outcome, muac, W_AZ, hg, Distancenew))#using only 
p_simple

########################################################################################################################################
fit.logistic.waz <-  glm(outcome ~ W_AZ,data=df, family=binomial(link = 'logit'))
pred.prob.logistic.waz <- predict(fit.logistic.waz,type='response')
plot(y=(pred.prob.logistic.waz),x=( df$W_AZ[-c(1:6)]))
temp <- df %>% select(W_AZ) %>% drop_na()
temp$pred.prob <- pred.prob.logistic.waz
temp$pred.prob <- round(temp$pred.prob, 3)
temp %>% filter(pred.prob == 0.037) %>% summarise(mean = mean(W_AZ), 
                                                  sd = sd(W_AZ))

fit.logistic.muac <- glm(outcome ~ muac,data=df, family=binomial(link = 'logit'))
pred.prob.logistic.muac <- predict(fit.logistic.muac,type='response')
temp <- df %>% select(muac) %>% drop_na()
temp$pred.prob <- pred.prob.logistic.muac
temp$pred.prob <- round(temp$pred.prob, 3)
temp %>% filter(pred.prob == 0.040) %>% summarise(mean = mean(muac), 
                                                  sd = sd(muac))

fit.logistic.edu <- glm(outcome ~ maternal_educationnew,data=df, family=binomial(link = 'logit'))
pred.prob.logistic.edu <- predict(fit.logistic.edu,type='response')

hg.df <- df %>% select(outcome,hg) %>% drop_na()
fit.logistic.hg <- glm(outcome ~ hg,data=hg.df, family=binomial(link = 'logit'))
pred.prob.logistic.hg <- predict(fit.logistic.hg,type='response')
temp <- df %>% select(hg) %>% drop_na()
temp$pred.prob <- pred.prob.logistic.hg
temp$pred.prob <- round(temp$pred.prob, 3)
temp %>% filter(pred.prob == 0.035) %>% summarise(mean = mean(hg), 
                                                  sd = sd(hg))
plot((temp$hg),temp$pred.prob)



########################################################################################################################################
#Decsion Tree learning 

df.small <- df %>% select(sex,W_AZ,temperature,spo2,Abnormal_BCS,time_since_last_hospitaliz,
                          hg,malaria,sibling_deaths,boiling,hivFinal,bednet,maternal_educationnew, Distancenew, hivMaternal,outcome)

K <- 10
cv.precsion <- cv.recall <- cv.fp_rat <- cv.auc <- cv.spec <- rep(0,K)
for(j in 1:K){
  
  #parameters for rpart
  loss.1 <- matrix(c(0,40,1,0), nrow=2)
  k <- 10
  parms <- list(split = 'information', loss =loss.1 )
  formula.fit <- as.formula("outcome ~ .")
  control <- rpart.control(minsplit = 40,minbucket = 10,cp = 0.022,maxdepth = 20)
  
  
  #fit tree
  cv_out <- cv_classTree(df = df.small, k = k, parms = parms, formula = formula.fit, control = control)
  
  #get preformance metrics
  mean_cv <- cv_out %>% summarise_all(., mean)
  cv.precsion[j] <- mean_cv$precsion; cv.recall[j] <- mean_cv$recall
  cv.fp_rat[j] <- mean_cv$fp_rat; cv.auc[j] <- mean_cv$auc; cv.spec[j] <- mean_cv$spec 
}

#plot preformance results
cv_df <- tibble(cv.precsion = cv.precsion,
                cv.recall = cv.recall,
                cv.fp_rat = cv.fp_rat,
                cv.auc = cv.auc, 
                cv.spec = cv.spec)
cv_df %>% 
  gather(.) %>% 
  ggplot(.,aes(x =key, y = value, fill= factor(key)))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_brewer(palette = 'Dark2')+
  facet_wrap( key ~ .,scales = "free")


cv_df %>% summarise_all(mean)

########################################################################################################################################
#use mlr package to find optimal parameters for DT
library(mlr)
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
  makeNumericParam("maxdepth", lower = 10, upper = 40)
)

ctrl.grid <- makeTuneControlGrid()
resamp.cv <- makeResampleDesc("CV", iters= 10)
set.seed(18)
tuned.mort <- tuneParams(lrn, task = mort.task, resampling = resamp.cv,
                         control = ctrl.grid, par.set = ps,measures = list(tpr, fnr))

opt.grid <- as.data.frame(tuned.mort$opt.path)

tuned.mort2 <- tuneParams(lrn, task = mort.task, resampling = resamp.cv,
                          control = ctrl.grid, par.set = ps,measures = list(fn))
########################################################################################################################################
#fit opt dt and print results
set.seed(100)
loss.1 <- matrix(c(0,40,1,0), nrow=2)
loss.2 <- matrix(c(0,45,1,0), nrow=2)
parms.1 <- list(split = 'information', loss =loss.1 )
parms.2 <- list(split = 'information', loss = loss.2)
formula.fit <- as.formula("outcome ~ .")
control <- rpart.control(minsplit = 40,minbucket = 10,cp = 0.022,maxdepth = 20)

tree.fit.1 <-  rpart(formula = formula.fit, data = df.small, method = 'class', 
                       parms = parms.1, control = control)
cp.min <- tree.fit.1$cptable[,1][which.min(tree.fit.1$cptable[,4])]
tree.pruned.1 <- prune.rpart(tree.fit.1, cp.min)
plot(tree.pruned.1);text(tree.pruned.1)

set.seed(102)
tree.fit.2 <- rpart(formula = formula.fit, data = df.small, method = 'class', 
                    parms = parms.2, control = control)
cp.min <- tree.fit.2$cptable[,1][which.min(tree.fit.2$cptable[,4])]
tree.fit.2 <- prune.rpart(tree.fit.2, cp.min)
plot(tree.fit.2);text(tree.fit.2)
########################################################################################################################################
#Jeff suggestion: Fit a rf then use importance rankings to decide what variables to include in DT. 
#after that find optimal parameters for DT and produce plot

rf.1 <- randomForest(outcome ~ ., data = df, ntree=1000, na.action = na.omit)
varImpPlot(rf.1)

rf.2 <- randomForest(outcome ~ ., data = df, ntree=2000, na.action = na.omit)
varImpPlot(rf.2)

rf.3 <- randomForest(outcome ~ ., data = df, ntree=5000, na.action = na.omit)
varImpPlot(rf.3)

rf.4 <- randomForest(outcome ~ ., data = df, ntree=10000, na.action = na.omit)
varImpPlot(rf.4)

#variable importance measure very different from LASSO
#interesting --> explore stability further
library(partykit)


rf.5 <- cforest(outcome ~ ., ntree = 1000,  data = df)
imp <- varimp(rf.5)



