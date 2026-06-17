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


#########################################################
######### predict stress in each plant species #########
#########################################################

# split the original object by species
ps_sp_l<-lapply(CSS_BacFun_rhizoplane_HighFilt_ps_l, function (x)
  phyloseq_sep_variable(x,variable = c("Sp_abb_name")))


#### Bacterial communities
#### MeJA VS control

set.seed(101)
Boruta_stress_bac_33sp<-mclapply(ps_sp_l$Bac, function (x) #parLapply wil do a parallel lapply on the defined cluster
  Boruta(Stress~.,   # classification you are trying to predict
         data = single_physeq_to_borutaInput(physeq_object = subset_samples(x, Stress %in% c("Control", "MeJA")),
                                             variable_to_be_classified = "Stress")[,-1], # removes first column, "sample" as a predictor variable
         doTrace=0,
         maxRuns = 1000,  #increase the maximum number of runs to decrease the number of tenttively important OTUs.
         ntree = 400000, # increase the number of trees to increase precision. decrease ntree/maxruns to reduce computational time.
         mc.cores = detectCores()))



save(Boruta_stress_bac_33sp, file = "./10_Machine_learning/r_objects/Boruta_MeJA_bac_33sp_400k.RData")


#### Fungal communities

set.seed(101)
Boruta_stress_fun_33sp<-mclapply(ps_sp_l$Fun, function (x) #parLapply wil do a parallel lapply on the defined cluster
  Boruta(Stress~.,
         data = single_physeq_to_borutaInput(physeq_object = subset_samples(x, Stress %in% c("Control", "MeJA")),
                                             variable_to_be_classified = "Stress")[,-1], # removes first column, "sample" as a predictor variable
         doTrace=0,
         maxRuns = 1000,  #increase the maximum number of runs to decrease the number of tenttively important OTUs.
         ntree = 400000, # increase the number of trees to increase precision. decrease ntree/maxruns to reduce computational time.
         mc.cores = detectCores()))


save(Boruta_stress_fun_33sp, file = "./10_Machine_learning/r_objects/Boruta_MeJA_fun_33sp_400k.RData")


#### Bacterial communities
#### SA VS control



set.seed(101)
Boruta_stress_bac_33sp<-mclapply(ps_sp_l$Bac, function (x) #parLapply wil do a parallel lapply on the defined cluster
  Boruta(Stress~.,   # classification you are trying to predict
         data = single_physeq_to_borutaInput(physeq_object = subset_samples(x, Stress %in% c("Control", "SA")),
                                             variable_to_be_classified = "Stress")[,-1], # removes first column, "sample" as a predictor variable
         doTrace=0,
         maxRuns = 1000,  #increase the maximum number of runs to decrease the number of tenttively important OTUs.
         ntree = 400000, # increase the number of trees to increase precision. decrease ntree/maxruns to reduce computational time.
         mc.cores = detectCores()))



save(Boruta_stress_bac_33sp, file = "./10_Machine_learning/r_objects/Boruta_SA_bac_33sp_400k.RData")


#### Fungal communities

set.seed(101)
Boruta_stress_fun_33sp<-mclapply(ps_sp_l$Fun, function (x) #parLapply wil do a parallel lapply on the defined cluster
  Boruta(Stress~.,
         data = single_physeq_to_borutaInput(physeq_object = subset_samples(x, Stress %in% c("Control", "SA")),
                                             variable_to_be_classified = "Stress")[,-1], # removes first column, "sample" as a predictor variable
         doTrace=0,
         maxRuns = 1000,  #increase the maximum number of runs to decrease the number of tenttively important OTUs.
         ntree = 400000, # increase the number of trees to increase precision. decrease ntree/maxruns to reduce computational time.
         mc.cores = detectCores()))


save(Boruta_stress_fun_33sp, file = "./10_Machine_learning/r_objects/Boruta_SA_fun_33sp_400k.RData")

