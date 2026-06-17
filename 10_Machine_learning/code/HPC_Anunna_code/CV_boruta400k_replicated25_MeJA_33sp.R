# load boruta_l_33sp, run fix_split_traintest_replicated function
# do this for 20k borutas and 100k borutas

# this will load all key libraries for the whole project
library(metagMisc)
library(dplyr)
library(tibble)
library(tidyverse)
library(Boruta)
library(phyloseq)
library(caret)
library(randomForest)
library(parallel) # use parallel computing to speed up RF calculations

# this will load the CSS 16S and ITS data as a list of two objects
load("./Data/phyloseq_objects/CSS_BacFun_rhizoplane_HighFilt_ps_l.RData")

# loads the custom functions used in this script
source("./10_Machine_learning/code/Machine_Learning_custom_functions.R")

######### split

# split the original object by species
ps_sp_l<-lapply(CSS_BacFun_rhizoplane_HighFilt_ps_l, function (x)
  phyloseq_sep_variable(x,variable = c("Sp_abb_name")))


########## load boruta objects

ps_sp_l$Bac
load(file = "./10_Machine_learning/r_objects/Boruta_MeJA_bac_33sp_400k.RData") #CHANGEd to 400k
Boruta_stress_bac_33sp

ps_sp_l$Fun
load(file = "./10_Machine_learning/r_objects/Boruta_MeJA_fun_33sp_400k.RData")
Boruta_stress_fun_33sp








#overwirte a few functions to make sure they are using mcmapply in the HPC

# write a function to split training and test data from a phyloseq object
train_and_test_spliter<-function(ps_object, variable_to_be_classified){
  
  # this function will separate train and test sets based on a phyloseq object and the variable to be predicted. it requires the function single_physeq_to_borutaInput()
  # ps_object = a phyloseq object
  # variable_to_be_classified =  a (quoted) metadata column that you want to predict
  # the output is a list of two objects: the first is the training set, the second is the test set
  
  # wrangle phyloseq data
  ps_data<-single_physeq_to_borutaInput(physeq_object = ps_object,
                                        variable_to_be_classified = variable_to_be_classified)
  
  # define training and test set. this can be ofptimized for repeated k-fold cross validation
  trainIndex<- createDataPartition(ps_data[,2],
                                   p = .70,
                                   list = FALSE,
                                   times = 1)
  # set train and test sets
  data_Train <- ps_data [ trainIndex,]
  data_Test  <- ps_data [-trainIndex,]
  
  output<-list(data_Train,data_Test)
  names(output)<-c("data_Train","data_Test")
  
  return(output)
  
}





# define a function to fix borta objects and put them into formula format
fixed_boruta_formula<-function(boruta_object){
  # this fucntion takes a boruta ofbect, fixes the inconclusive tas into importnat o unimportnat, and then generates a formula
  # the input is a boruta object
  # the output is a boruta formula to be fed to caret::train
  # NOTE: boruta objects with zero imporntat features may crash!
  
  fixed_boruta<-TentativeRoughFix(boruta_object)
  boruta_imp_ASV<-getSelectedAttributes(fixed_boruta)
  print("number of importnat ASVs. Warning: if zero, formula will crash!")
  print(length(boruta_imp_ASV)%>%unlist()%>%sort())
  formula_boruta<-getConfirmedFormula(fixed_boruta)
  
  return(formula_boruta)
}



# defines custom fucntion to fix, split train and test data
fix_split_train_test_replicated<-function (boruta_output_l, ps_object_l, variable_to_be_classified){
  
  
  # fix boruta in a formula to be evaluated with caret
  boruta_formula_bac_l<-lapply(boruta_output_l, function(x) fixed_boruta_formula(x))
  
  # split train adn test dataset
  train_test_l<-mclapply(ps_object_l, function (x)
    train_and_test_spliter(ps_object = x,
                           variable_to_be_classified = variable_to_be_classified),
    mc.cores = detectCores())
  
  
  
  # train model
  boruta_feature_rf_repeatedcv<-mcmapply(function (x,z) {
    
    train.control <- caret::trainControl(method = "repeatedcv", # set trainig/data split controls for the train function
                                         number = 5,
                                         repeats = 50,
                                         allowParallel = TRUE,
                                         returnData = FALSE)
    
    model_borutized <- caret::train(form = z, # bruta formula
                                    data = x[[1]], # training data ; first element of train_and_test_spliter()
                                    method = "rf", #execute training based on RF
                                    trControl = train.control, # defined in trainControl() above
                                    ntree=8000)
    
    
    
    return(model_borutized)
  },
  x = train_test_l,
  z = boruta_formula_bac_l,
  mc.cores = 4, # trying to limit this to avoid OMM-killing
  SIMPLIFY = FALSE)
  
  
  
  #test model
  confusion_matrix_output<-mcmapply(function(x,y){
    prediction<-stats::predict(object = x, newdata = y[[2]])
    confusion_output<-confusionMatrix(data = prediction, reference = y[[2]][,2])
    return(confusion_output)
  },
  x = boruta_feature_rf_repeatedcv,
  y = train_test_l,
  mc.cores = 4,
  SIMPLIFY = FALSE)
  
  
  
  output<-list("trained_model_rf_repeatedcv" = boruta_feature_rf_repeatedcv,
               "confusion_matrix_output" = confusion_matrix_output)
  
  return(output)
  
  
}






#train bacteria

t0<-Sys.time()
set.seed(23456)
replicated_model_400k_bac_MeJA<-replicate(25, 
                                 fix_split_train_test_replicated(boruta_output_l = Boruta_stress_bac_33sp,
                                                                 ps_object_l = lapply(ps_sp_l$Bac, function (x)
                                                                   subset_samples(x,
                                                                                  Stress %in% c("Control",
                                                                                                "MeJA"))),
                                                                 variable_to_be_classified ="Stress"))


t1<-Sys.time()
t0-t1



# now train fungi

t0<-Sys.time()
set.seed(23456)
replicated_model_400k_fun_MeJA<-replicate(25, 
                                 fix_split_train_test_replicated(boruta_output_l = Boruta_stress_fun_33sp[-31], # remove maize as it has zero important ASVs
                                                                 ps_object_l = lapply(ps_sp_l$Fun[-31], function (x)
                                                                   subset_samples(x,
                                                                                  Stress %in% c("Control",
                                                                                                "MeJA"))),
                                                                 variable_to_be_classified ="Stress"))


t1<-Sys.time()
t0-t1



CV_25replicated_MeJA_400k_l_l_33sp<-list("replicated_model_400k_bac_MeJA" = replicated_model_400k_bac_MeJA,
                                         "replicated_model_400k_fun_MeJA" = replicated_model_400k_fun_MeJA)

save(CV_25replicated_MeJA_400k_l_l_33sp, file = "./10_Machine_learning/r_objects/CV_25replicated_MeJA_400k_l_l_33sp.RData")
