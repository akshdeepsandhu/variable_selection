#funcs_analysis
#functions that will be used for anlaysis 
#one file with functions to keep Rmd clean 

pre_processing <- function(df, del.vars,factor.vars){
  #Function that will preprocess the data for model fitting
  
  #change var types and delete variables
  df <- df %>% 
    dselect(-del.vars) %>% 
    mutate_at(vars(factor.vars), funs(factor))

  #centre and scale 
  
  
  
  
  
  return(df)
  
}
