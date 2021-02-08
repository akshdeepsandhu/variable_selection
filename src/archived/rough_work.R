#Rough work
mort0 <- 1325
mort1 <- 70
#######################################################################################################################
#Funcs.R

cutOff_acc <- function(cut_off,var,outcome){
  #Function that fits logistic regression, 
  #then looks at predicted probabilites and determines clas, 
  #based on cut_off parameter
  #Returns tp and fp (true and false postive numbers)
  
  
  fit.logit <- glm(outcome ~ var, family = binomial(link='logit'))
  fitted <- fitted(fit.logit)
  cut_off_df <- tibble(outcome = outcome, prob = fitted, var = var)
  cut_off_df <- cut_off_df %>% mutate(pred_class = if_else(prob > cut_off, 1,0))
  #get number of true postives (tp) and false postives (fp)
  cut_off_df <- cut_off_df %>% mutate(tp = if_else(outcome == 1 & pred_class == 1, 1,0 ))
  cut_off_df <- cut_off_df %>% mutate(fp = if_else(outcome == 0 & pred_class == 1, 1, 0))
  
  tp <- sum(cut_off_df$tp)
  fp <- sum(cut_off_df$fp)
  
  
  
  
  return_vec <- c(tp,fp)
  
  return(return_vec)
  
}

plot_cutOff_results <- function(outcome,var,grid, class0, class1){
  #recall -> class is correctly recognised (tp / sum(p))
  #precision -> postive actually postive (tp/ tp + fp)
  #High recall, low precsion -> lots of postives correctly recongnized (low fn) but lots of fp
  #low recall, high precsion -> high fn, but low fp
  
  output <- sapply(grid, cutOff_acc, var = var, outcome = outcome)
  output.df <-  tibble(recall = output[1,]/class1, precsion = output[1,]/(output[1,]+output[2,]), fp_pct =output[2,]/class0 ,cut_off = grid)
  #plot precsion and recall
  p1 <- output.df%>% 
    ggplot(.,aes(x=cut_off))+
    #  geom_line()+
    geom_line(aes(y=recall,colour='recall'))+
    geom_line(aes(y=precsion,colour='precsion'))+
    geom_line(aes(y=fp_pct,colour='fp_pct'))+
    theme_bw()+
    ylab('rate')
  #    scale_y_continuous(limits=c(0,1))+
  #    scale_x_continuous(limits=c(0,1))
  
  
  #print and return plot and df
  print(p1)
  return(list(p1, output.df))
  
}


#######################################################################################################################
#Funcitons above fit logistic regression and plot precsion-recall curves 
#explore precsion recall curves for all explanatory vairables in df
#recall -> class is correctly recognised (tp / sum(p))
#precision -> postive actually postive (tp/ tp + fp)
#High recall, low precsion -> lots of postives correctly recongnized (low fn) but lots of fp
#low recall, high precsion -> high fn, but low fp
plot_cutOff_results(outcome = df$final_outcome_death,var = df$sex,grid = seq(0,0.06,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$AgeAtAdmiImp,grid = seq(0,0.2,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$muac,grid = seq(0,0.2,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$W_AZ,grid = seq(0,0.3,0.005),class0 = mort0,class1 = mort1)#decent candidate
plot_cutOff_results(outcome = df$final_outcome_death,var = df$W_LZ,grid = seq(0,0.15,0.005),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$BMI_AZ,grid = seq(0,0.15,0.005),class0 = mort0,class1 = mort1)#not bad
plot_cutOff_results(outcome = df$final_outcome_death,var = df$hr,grid = seq(0,0.075,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$rr,grid = seq(0,0.075,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$sbp,grid = seq(0,0.075,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$dbp,grid = seq(0,0.075,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$temperature,grid = seq(0,0.1,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$spo2,grid = seq(0,0.1,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$Abnormal_BCS,grid = seq(0,0.2,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$time_since_last_hospitaliz,grid = seq(0,0.075,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$malaria,grid = seq(0,0.075,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$number_of_siblings,grid = seq(0,0.075,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$sibling_deaths,grid = seq(0,0.075,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$maternal_age_years,grid = seq(0,0.075,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$boiling,grid = seq(0,0.1,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$hivFinal,grid = seq(0,0.2,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$bednet,grid = seq(0,0.075,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$bednet,grid = seq(0,0.1,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$maternal_educationnew,grid = seq(0,0.1,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$Distancenew,grid = seq(0,0.1,0.001),class0 = mort0,class1 = mort1)#decent candidate
plot_cutOff_results(outcome = df$final_outcome_death,var = df$hivMaternal,grid = seq(0,0.2,0.001),class0 = mort0,class1 = mort1)
plot_cutOff_results(outcome = df$final_outcome_death,var = df$water_source_drinking,grid = seq(0,0.1,0.001),class0 = mort0,class1 = mort1)
#all variables have fairly similar curves, starting with high recall and low precsion. Precsion usally stays very low, 
#so it might make sense to ignore it and try to find a cutoff that has a high fp rate but good recall. 
#need to clarify with Jeff 