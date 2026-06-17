



#prepare phyloseq object to be an input in Boruta
# the first input argument "list_physeq_object" is a list of phyloseq objects
# se second input argumnet "variable_to_be_classified" is a quoted variable, present in your phyloseq metadata, which you want to predict in random forest

physeq_to_borutaInput<-function (list_physeq_object, variable_to_be_classified){
  # boruta expects a transposed OTU table with variable_to_be_classified as a added variable
  # the output is a list of df ready to be used as input to boruta  
  #transpose phtseq otu table  
  otu_cells_wide_list <-lapply(list_physeq_object, function(x) #transpose the feature table...
    base::as.data.frame(t(otu_table(x)))%>%
      rownames_to_column(var = "sample"))
  
  # extract sample data
  metadata_list <-lapply(list_physeq_object, function(x)
    as(sample_data(x),"data.frame")%>%
      rownames_to_column(var = "sample"))
  
  #add the variable classification you want to predict with the random forest
  boruta_dataset<-mapply(function (x,y)
    base::merge(dplyr::select(x,sample, variable_to_be_classified),
                y,
                by = "sample",
                all.y = TRUE),
    x = metadata_list,
    y = otu_cells_wide_list,
    SIMPLIFY = FALSE)
  
  #make sure your variable to be classified is a factor, or boruta won't run
  output<-lapply(boruta_dataset, function(x) {
    x[,2]<-as.factor(x[,2]) # saves the second column, your variable_to_be_classified, as a factor
    return(x)
  })
  gc()
  return(output)
}












single_physeq_to_borutaInput<-function (physeq_object, variable_to_be_classified){
  # boruta expects a transposed OTU table with variable_to_be_classified as a added variable
  # the output is a list of df ready to be used as input to boruta  
  #transpose phtseq otu table  
  otu_cells_wide_list <- #transpose the feature table...
    base::as.data.frame(t(otu_table(physeq_object)))%>%
      rownames_to_column(var = "sample")
  
  # extract sample data
  metadata_list <-
    as(sample_data(physeq_object),"data.frame")%>%
      rownames_to_column(var = "sample")
  
  #add the variable classification you want to predict with the random forest
  boruta_dataset<-
    base::merge(dplyr::select(metadata_list,sample, variable_to_be_classified),
                otu_cells_wide_list,
                by = "sample",
                all.y = TRUE)
  
  #make sure your variable to be classified is a factor, or boruta won't run
  
    boruta_dataset[,2]<-as.factor(boruta_dataset[,2]) # saves the second column, your variable_to_be_classified, as a factor
   
    output<-boruta_dataset
    
  gc()
  return(output)
}






# define a funciton that obtains the samples closes to centroids
# ultimatelty, this function won't be useful for our random forests, but might be handy later
get_nn_centroid_samples<- function(ordinate_output, nn, factor_column_n, ps_object){
  
  # this function will  define centroids, calculate distances, and extract sample names of the nn samples closes to the centroid
  # ordinate_output = output form the phyloseq::ordinate() function. this was only tested for bray-curtis distances
  # nn = number of samples closest to centroid you want to extract (such as 3 closest samples) 
  # factor_column_n = number of the column with relevant emtadata you want to extract (such as sp_and_stress, column 46 of our metadata)
  # ps_object = a phyloseq objectec, used to calculate the ordinations and that contains the emtadata columns
  
  
  # save the points of the ordination as a df
  NMDS_coordinates<-as.data.frame(ordinate_output$points)
  
  #merge these coordinates with the rest of the metadata
  NMDS_coordinates<-merge(as.data.frame(ps_object@sam_data),
                          NMDS_coordinates,
                          by = 0)
  
  
  # split coordinates by plant sp and stress
  split_points<-split(NMDS_coordinates, f = NMDS_coordinates[,factor_column_n+1]) # this +1 refers to changes in column names due to the merging that moved row names to a column
  
  # calculate centroid of each sp + stress conditions
  split_points<-lapply(split_points, function(x)dplyr::mutate(x, mean_MDS1 = mean(MDS1), mean_MDS2 = mean(MDS2)))
  
  # remove unecessary emtadata
  split_points<-lapply(split_points, function(x) {
    rownames(x)<-x[,1]
    #x<-x [,47:50]
    
    x<-x[, c("MDS1", "MDS2", "mean_MDS1", "mean_MDS2")]
    return(x)
  })
  
  
  
  # this adds the euclediandistance to the control centroid, species by species and stress by stress
  NMDS_centroids_l<-lapply(split_points, function (x){ # for each species+stress of the list...
    hypotenuse<-apply(x,1, function(y) # perform the following function for every row (only on numeric cols of the input df)
      sqrt(sum((x[1,3]-y[1])^2)+sum((x[1,4]-y[2])^2))) 
    x$distance_to_centroid<-hypotenuse
    return(x)
    
    
    # x[1,3] = NMDS1_centroid; 
    # y[1] = NMDS1_row; 
    # x[1,4] NMDS2_centroid;; 
    # y[2] = NMDS2_row; 
    # elevate on power 2 and then take a square root to destroy the signal
    # these two values will be the sides of your triangule
    # elevate both sides of the triagule to power 2 and then take a square root to find the hypotenuse
    # the hypotenuse is the distance between the centroids
  })
  
  # sort sampels according distance to centroid
  NMDS_centroids_l<-lapply(NMDS_centroids_l, function(x) x[order(x$distance_to_centroid, decreasing=FALSE),])
  
  # get the 3 samples that are closest to centroid
  NMDS_centroids_l<-lapply(NMDS_centroids_l, function(x) rownames(x[1:nn,]))
  
  #define output
  output<-NMDS_centroids_l
  
  return(output)
  
}
















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


# put fixing, spliting, training and testing all in a single function
fix_split_train_test<-function (boruta_output_l, ps_object_l, variable_to_be_classified){
  # this fucntion will fix tentative features in a list of boruta objects
  # then it will split a list of phyloseq objects into training and test sets
  # then it will train the list of models
  # then it will test the lsit of models
  # then it returns a list of confusion matrixes (one for each model)
    # boruta_output_l = a list of boruta objects
    # ps_object_l = a list of phyloseq objects
    # variable_to_be_classified = the metadata varaible you are trying to predict (must be quoted, like "Stress")
  
  # NOTE ON SETTING SEED: 
  # if set.seed(3456) is kept INSIDE the function, results will be identical to another run where set.seed(3456) is kept inside the function
  # if set.seed(3456) is kept OUTSIDE the function, silencing the sed.seed inside of it, results will be identical to another run where set.seed(3456) is kept OUTSIDE the function
  # a run with set.seet(3456) INSIDE the function will be different from a run with set.seet(3456) OUTSIDE the function
  # when the function is replicated, set.seet set OUTSIDE the function will work well (different traint/test data, predictions on the replicated set; kept  consistent under the same seed for the set); time taken to calculate models will differ
  
  # fix boruta in a formula to be evaluated with caret
  boruta_formula_bac_l<-lapply(boruta_output_l, function(x) fixed_boruta_formula(x))
  
  # split train adn test dataset
  #set.seed(3456)
  train_test_l<-lapply(ps_object_l, function (x)
    train_and_test_spliter(ps_object = x, 
                           variable_to_be_classified = variable_to_be_classified))
  
  
  
  # train model
  #set.seed(3456)
  boruta_feature_rf_repeatedcv<-mapply(function (x,z) {
    
    train.control <- caret::trainControl(method = "repeatedcv", # set trainig/data split controls for the train function
                                         number = 5,
                                         repeats = 100,
                                         allowParallel = TRUE)
    
    model_borutized <- caret::train(form = z, # bruta formula
                                    data = x[[1]], # training data ; first element of train_and_test_spliter()
                                    method = "rf", #execute training based on RF
                                    trControl = train.control, # defined in trainControl() above
                                    ntree=10000)
    
    
    
    return(model_borutized)
  },
  x = train_test_l,
  z = boruta_formula_bac_l,
  SIMPLIFY = FALSE)
  
  
  
  #test model
 # set.seed(3456)
  confusion_matrix_output<-mapply(function(x,y){
    prediction<-stats::predict(object = x, newdata = y[[2]]) 
    confusion_output<-confusionMatrix(data = prediction, reference = y[[2]][,2])
    return(confusion_output)
  },
  x = boruta_feature_rf_repeatedcv,
  y = train_test_l,
  SIMPLIFY = FALSE)
  
  output<-list("trained_model_rf_repeatedcv" = boruta_feature_rf_repeatedcv,
               "confusion_matrix_output" = confusion_matrix_output)
  
  return(output)
  
  
}






# put fixing, spliting, training and testing all in a single function, to be replicated multiple times
fix_split_train_test_replicated<-function (boruta_output_l, ps_object_l, variable_to_be_classified){
  # this fucntion will fix tentative features in a list of boruta objects
  # then it will split a list of phyloseq objects into training and test sets
  # then it will train the list of models
  # then it will test the lsit of models
  # then it returns a list of confusion matrixes (one for each model)
  # boruta_output_l = a list of boruta objects
  # ps_object_l = a list of phyloseq objects
  # variable_to_be_classified = the metadata varaible you are trying to predict (must be quoted, like "Stress")
  
  # NOTE ON SETTING SEED: 
  # if set.seed(3456) is kept INSIDE the function, results will be identical to another run where set.seed(3456) is kept inside the function
  # if set.seed(3456) is kept OUTSIDE the function, silencing the sed.seed inside of it, results will be identical to another run where set.seed(3456) is kept OUTSIDE the function
  # a run with set.seet(3456) INSIDE the function will be different from a run with set.seet(3456) OUTSIDE the function
  # when the function is replicated, set.seet set OUTSIDE the function will work well (different traint/test data, predictions on the replicated set; kept  consistent under the same seed for the set); time taken to calculate models will differ
  
  # fix boruta in a formula to be evaluated with caret
  boruta_formula_bac_l<-lapply(boruta_output_l, function(x) fixed_boruta_formula(x))
  
  # split train adn test dataset
  #set.seed(3456)
  train_test_l<-lapply(ps_object_l, function (x)
    train_and_test_spliter(ps_object = x, 
                           variable_to_be_classified = variable_to_be_classified))
  
  
  
  # train model
  #set.seed(3456)
  boruta_feature_rf_repeatedcv<-mapply(function (x,z) {
    
    train.control <- caret::trainControl(method = "repeatedcv", # set trainig/data split controls for the train function
                                         number = 5,
                                         repeats = 50,
                                         allowParallel = TRUE)
    
    model_borutized <- caret::train(form = z, # bruta formula
                                    data = x[[1]], # training data ; first element of train_and_test_spliter()
                                    method = "rf", #execute training based on RF
                                    trControl = train.control, # defined in trainControl() above
                                    ntree=8000)
    
    
    
    return(model_borutized)
  },
  x = train_test_l,
  z = boruta_formula_bac_l,
  SIMPLIFY = FALSE)
  
  
  
  #test model
  # set.seed(3456)
  confusion_matrix_output<-mapply(function(x,y){
    prediction<-stats::predict(object = x, newdata = y[[2]]) 
    confusion_output<-confusionMatrix(data = prediction, reference = y[[2]][,2])
    return(confusion_output)
  },
  x = boruta_feature_rf_repeatedcv,
  y = train_test_l,
  SIMPLIFY = FALSE)
  
  output<-list("trained_model_rf_repeatedcv" = boruta_feature_rf_repeatedcv,
               "confusion_matrix_output" = confusion_matrix_output)
  
  return(output)
  
  
}










# define function to extract some key model metrics after CV precision testing
extract_confusionmatrix_metrics<-function(CV_output_l, Imp_ASV_stats_l){
  # this fucntion will extract key RF model metris from lists of RF models
  # CV_output_l = list of model testing objects from fix_split_train_test() custom function
  # Imp_ASV_stats_l = list of important ASV stats from a fixed boruta object
  # the output is a df with metrics for each model

  # accurayc, kappa, AccuracyPValue  
  
  cv_metrics<-map(CV_output_l,3)
  
  accuracy<-unlist(map(cv_metrics,1)) #accuracy
  
  kappa<-unlist(map(cv_metrics,2)) #kappa
  AccuracyLower<-unlist(map(cv_metrics,3)) #AccuracyLower  
  AccuracyUpper<-unlist(map(cv_metrics,4)) #AccuracyUpper  
  AccuracyPValue<-unlist(map(cv_metrics,6)) #AccuracyPValue 
  
  
  
  model_metrics<-as.data.frame(accuracy)
  model_metrics$kappa<-kappa     
  model_metrics$AccuracyLower<-AccuracyLower  
  model_metrics$AccuracyUpper<-AccuracyUpper  
  model_metrics$AccuracyPValue<-AccuracyPValue
  
  # n important stress predictor ASVs  
  model_metrics$n_imp_ASVs<-lapply(Imp_ASV_stats_l, function (x) length(x[,1]))%>%unlist
  
  # mean average imporance of ASVs in model
  model_metrics$mean_ASV_imp<-lapply(Imp_ASV_stats_l, function (x) mean(x[,2]))%>%unlist
  
  # mean median imporance of ASVs in model
  model_metrics$median_ASV_imp<-lapply(Imp_ASV_stats_l, function (x) mean(x[,3]))%>%unlist
  
  # sd imporance of ASVs in model
  model_metrics$sd_ASV_imp<-lapply(Imp_ASV_stats_l, function (x) sd(x[,2]))%>%unlist
  
  # men normHits of ASVs in model
  model_metrics$mean_normHits_ASV_imp<-lapply(Imp_ASV_stats_l, function (x) mean(x[,6]))%>%unlist
  
  # sd normHits of ASVs in model
  model_metrics$sd_normHits_ASV_imp<-lapply(Imp_ASV_stats_l, function (x) sd(x[,6]))%>%unlist
  
  return(model_metrics)
}
