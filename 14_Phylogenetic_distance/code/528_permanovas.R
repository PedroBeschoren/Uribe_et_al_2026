#load libraries
library("vegan")
library("phyloseq")
library("dplyr")
library("parallel")


####### Bac

# load output of pipeline
load(file ="./14_Phylogenetic_distance/r_objects/bac_alltaxa_phylosig_output_PCoA_nonsig.RData")


# testing sp*stress per taxonomic group is a priority. this function will perform that
permanova_per_taxonomic_group_spXstress<-function(ps_obj_l){
  # ps_obj_l = a list of phyloseq objects, split by taxonomic group
  
  
  
  ps_obj_l<- lapply(ps_obj_l, function(x){
    
    x@sam_data$Sp_Lineage_Walden_extended<-
      if_else(condition = x@sam_data$Sp_abb_name %in% c("Ia", "Lm", "Co", "Bi"), 
              true ="lineage_II_extended",
              false =  x@sam_data$Sp_Lineage_Walden)
    
    return(x)
    
  })
  
  
  #calculate permanova at individual stress conditions
  set.seed(303848)
  permanova_l <- mclapply(ps_obj_l, function (x)
    adonis2(formula = phyloseq::distance(t(otu_table(x)), method="bray") 
            ~ Sp_full_name + Stress + Sp_full_name*Stress + Block,
            data = as(sample_data(x),"data.frame"), # changing with as.data.frame is insufficient
            permutations = 999, by = "terms"),  
    mc.cores = detectCores())
  
  
  beta_disper_l <- mclapply(ps_obj_l, function (x){
    
    beta_disper_sp <-  betadisper(phyloseq::distance(t(otu_table(x)), method = "bray"), 
                                  sample_data(x)$Sp_abb_name)
    beta_disper_stress <-  betadisper(phyloseq::distance(t(otu_table(x)), method = "bray"),
                                      sample_data(x)$Stress)
    beta_disper_block <-  betadisper(phyloseq::distance(t(otu_table(x)), method = "bray"),
                                     sample_data(x)$Block)
    
    output<-list("beta_disper_sp" = beta_disper_sp,
                 "beta_disper_stress" = beta_disper_stress,
                 "beta_disper_block" = beta_disper_block)
    
    return(output)     
  },  
  mc.cores = detectCores())
  
  
  
  output<-list("permanova_l" = permanova_l,
               "beta_disper_l" = beta_disper_l)
  
  return(output)
  
}

#run permanvoas
permanovas_spXstress_528<-
  permanova_per_taxonomic_group_spXstress(ps_obj_l = bac_alltaxa_phylosig_output_PCoA_nonsig$ordination_coordinates_l)


# save permanovas
save(permanovas_spXstress_528, file ="./14_Phylogenetic_distance/r_objects/permanovas_spXstress_528.RData")




# testing lineage per taxonomic group 
permanova_per_taxonomic_group_lineageXstress<-function(ps_obj_l){
  # ps_obj_l = a list of phyloseq objects, split by taxonomic group
  
  
  
  ps_obj_l<- lapply(ps_obj_l, function(x){
    
    x@sam_data$Sp_Lineage_Walden_extended<-
      if_else(condition = x@sam_data$Sp_abb_name %in% c("Ia", "Lm", "Co", "Bi"), 
              true ="lineage_II_extended",
              false =  x@sam_data$Sp_Lineage_Walden)
    
    return(x)
    
  })
  
  
  #calculate permanova at individual stress conditions
  set.seed(303848)
  permanova_l <- mclapply(ps_obj_l, function (x)
    adonis2(formula = phyloseq::distance(t(otu_table(x)), method="bray") 
            ~ Sp_Lineage_Walden_extended + Stress + Sp_Lineage_Walden_extended*Stress + Block,
            data = as(sample_data(x),"data.frame"), # changing with as.data.frame is insufficient
            permutations = 999, by = "terms"),  
    mc.cores = detectCores())
  
  
  beta_disper_l <- mclapply(ps_obj_l, function (x){
    
    beta_disper_Sp_Lineage_Walden_extended <-  betadisper(phyloseq::distance(t(otu_table(x)), method = "bray"), 
                                  sample_data(x)$Sp_Lineage_Walden_extended)
    beta_disper_Sp_Lineage_Walden <-  betadisper(phyloseq::distance(t(otu_table(x)), method = "bray"),
                                      sample_data(x)$Sp_Lineage_Walden)
   
    output<-list("beta_disper_Sp_Lineage_Walden_extended" = beta_disper_Sp_Lineage_Walden_extended,
                 "beta_disper_Sp_Lineage_Walden" = beta_disper_Sp_Lineage_Walden)
    
    return(output)     
  },  
  mc.cores = detectCores())
  
  
  
  output<-list("permanova_l" = permanova_l,
               "beta_disper_l" = beta_disper_l)
  
  return(output)
  
}

#run permanvoas
permanovas_lineageXstress_528<-
  permanova_per_taxonomic_group_lineageXstress(ps_obj_l = bac_alltaxa_phylosig_output_PCoA_nonsig$ordination_coordinates_l[20:21])


# save permanovas
save(permanovas_lineageXstress_528, file ="./14_Phylogenetic_distance/r_objects/permanovas_lineageXstress_528.RData")


####### Fun

# load output of pipeline
load(file ="./14_Phylogenetic_distance/r_objects/fun_alltaxa_phylosig_PCoA_nonsig.RData")


#run permanvoas
permanovas_spXstress_70<-
  permanova_per_taxonomic_group_spXstress(ps_obj_l = fun_alltaxa_phylosig_PCoA_nonsig$ordination_coordinates_l)

# save permanovas
save(permanovas_spXstress_70, file ="./14_Phylogenetic_distance/r_objects/permanovas_spXstress_70.RData")




#run permanvoas
permanovas_lineageXstress_70<-
  permanova_per_taxonomic_group_lineageXstress(ps_obj_l = fun_alltaxa_phylosig_PCoA_nonsig$ordination_coordinates_l)

# save permanovas
save(permanovas_lineageXstress_70, file ="./14_Phylogenetic_distance/r_objects/permanovas_lineageXstress_70.RData")


