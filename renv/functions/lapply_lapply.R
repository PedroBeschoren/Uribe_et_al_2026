# This function allows to run a custom function in a nested list. It is useful when working with our phyloseq object, as it contains a list of 2 (bacteria and fungia), each of which is made by a list of 37 samples (plant species + soil samples)
## see na example on chunk X.1.6 lists of lists


lapply_lapply<-function(nested_list,
                        function_for_single_object,
                        ...){
  
  # there are two inputs for this function
  # nested list = is a list of lists. usually a phyloseq objects with $bac and $fun. inside $bac you will find 37 plant species (a list of 37 things)
  # custom function  = a function that requires a single object as an argument, such as the custom function add_diversity_to_physeq_object
  # ... = arguments to be passed to the function_for_single_object, allowing you to adjust variables and paramenters
  
  # the output is the regular output for the function that requires a single object as an argument, but in indexed list format (list of 2, each part of the two has a list of 37)
  output<- lapply(nested_list, function(list_of_two)
    lapply(list_of_two, function (list_of_37)
      function_for_single_object(list_of_37, ...)))
  return(output)
}
