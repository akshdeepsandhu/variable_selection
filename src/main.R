##################
#Main func 
source('src/utils.R')
df <- data.create_matrix(data.clean_data(path = 'data/SD_for_Ash.csv'))[,-1]
n <- c(500,1000,2000,2500,3000,3500,4000)
samples.list <- data.generate_samples(df,n)
train.list <- samples.list[[1]]
test.df <- samples.list[[2]]


#variable selection frequency matrix
frequency.matrix <- data.frame(matrix(ncol = ncol(df) + 3, nrow = 0))
col_names <- colnames(df)[-c(1,20)]
colnames(frequency.matrix) <- c(col_names,"method", "N", "Spec4Sens", "AUC.TEST", "ROC.CV")

#set up clusters for parallel processing
cl <- makeCluster(detectCores())
clusterEvalQ(cl, {library(caret)
                library(broom)})

#fit rf w/ rfe model
model.list <- parLapply(cl,train.list, randomforest.fit_model)
metrics.list <- parLapply(cl,model.list, randomforest.model_metrics,test.df = test.df)
#add results to freq.matrix
frequency.matrix <- add_to_matrix(model.list,metrics.list,col_names,n,frequency.matrix,"RF")

#fit stepwise model
model.list <- parLapply(cl,train.list, stepwise.fit_model)
metrics.list <- parLapply(cl,model.list, stepwise.model_metrics,test.df = test.df)
#add results to freq.matrix
frequency.matrix <- add_to_matrix(model.list,metrics.list,col_names,n,frequency.matrix,"stepwiseAIC")

#fit univarate.pval model
model.list <- parLapply(cl,train.list, univariate.fit_model)
metrics.list <- parLapply(cl,model.list, univariate.model_metrics,test.df = test.df)
#add results to freq.matrix
frequency.matrix <- add_to_matrix(model.list,metrics.list,col_names,n,frequency.matrix,"univarate.pval")

#fit LASSO model
model.list <- parLapply(cl,train.list, lasso.fit_model)
metrics.list <- parLapply(cl,model.list, lasso.model_metrics,test.df = test.df)
#add results to freq.matrix
frequency.matrix <- add_to_matrix(model.list,metrics.list,col_names,n,frequency.matrix,"LASSO")

#fit RFE model
model.list <- parLapply(cl,train.list, rfe.fit_model)
metrics.list <- parLapply(cl,model.list, rfe.model_metrics,test.df = test.df)
#add results to freq.matrix
frequency.matrix <- add_to_matrix(model.list,metrics.list,col_names,n,frequency.matrix,"RFE")

#fit EN model
model.list <- parLapply(cl,train.list, en.fit_model)
metrics.list <- parLapply(cl,model.list, en.model_metrics,test.df = test.df)
#add results to freq.matrix
frequency.matrix <- add_to_matrix(model.list,metrics.list,col_names,n,frequency.matrix,"EN")

#fit adaLASSO model
model.list <- parLapply(cl,train.list, adaLasso.fit_model)
metrics.list <- parLapply(cl,model.list, adaLasso.model_metrics,test.df = test.df)
#add results to freq.matrix
frequency.matrix <- add_to_matrix(model.list,metrics.list,col_names,n,frequency.matrix,"adaLasso")




