#load libraries
library("vegan")
library("phyloseq")
library("dplyr")
library("parallel")
library("metagMisc")

####### Bac

# load output of pipeline, Bac
load(file ="./14_Phylogenetic_distance/r_objects/bac_alltaxa_phylosig_output_PCoA_nonsig.RData")


# load output of pipeline, Fun
load(file ="./14_Phylogenetic_distance/r_objects/fun_alltaxa_phylosig_PCoA_nonsig.RData")


#this function will extract the varaince each PCoA explained, and also plot 500 PCoAs
explained_PCoA_Axis_per_taxonomic_group<-function(ps_obj_l){
  # ps_obj_l = a list of phyloseq objects, split by taxonomic group
  
  # add extended lineage as metadata
  ps_obj_l<- lapply(ps_obj_l, function(x){
    
    x@sam_data$Sp_Lineage_Walden_extended<-
      if_else(condition = x@sam_data$Sp_abb_name %in% c("Ia", "Lm", "Co", "Bi"), 
              true ="lineage_II_extended",
              false =  x@sam_data$Sp_Lineage_Walden)
    
    return(x)
    
  })
  
  
  

  PCoA_l<- mclapply(ps_obj_l, function(x)
      ordinate(physeq = x,
               method="PCoA",
               distance="bray",
               autotransform=TRUE),
       mc.cores = detectCores())
    

  ploted_ordination<- mcmapply(function (x,y)
      plot_ordination(physeq = x,
                      ordination = y, 
                      color = "Sp_Lineage_Walden_extended",
                      shape = "Stress")+
      scale_shape_manual(values = c(19,23,3)),
      x = ps_obj_l,
      y = PCoA_l,
      SIMPLIFY = FALSE,
      mc.cores = detectCores())
    
  explained_variance<-lapply(ploted_ordination, function(z)
      c(z$labels$x, z$labels$y))

  output<-list("PCoA_l" = PCoA_l,
               "ploted_ordination" = ploted_ordination,
               "explained_variance" = explained_variance)
  
  return(output)
  
  
}


# run the PCoA axis pct and plot, Bac
PCoA_528_plots_axis_pct<-
  explained_PCoA_Axis_per_taxonomic_group(ps_obj_l = bac_alltaxa_phylosig_output_PCoA_nonsig$ordination_coordinates_l)

# save 528 PCoA plots & ordiantions, Bac
save(PCoA_528_plots_axis_pct, file = "./14_Phylogenetic_distance/r_objects/PCoA_528_plots_axis_pct.RData")

# run the PCoA axis pct and plot, Fun
PCoA_Fun_528_plots_axis_pct<-
  explained_PCoA_Axis_per_taxonomic_group(ps_obj_l = fun_alltaxa_phylosig_PCoA_nonsig$ordination_coordinates_l)

# save 528 PCoA plots & ordiantions, Fun
save(PCoA_Fun_528_plots_axis_pct, file = "./14_Phylogenetic_distance/r_objects/PCoA_Fun_528_plots_axis_pct.RData")




#split by stress, re-run function, Bac
Control_ps_l<-
  lapply(bac_alltaxa_phylosig_output_PCoA_nonsig$ordination_coordinates_l, function (x)
    subset_samples(x, Stress =="Control"))

MeJA_ps_l<-
  lapply(bac_alltaxa_phylosig_output_PCoA_nonsig$ordination_coordinates_l, function (x)
    subset_samples(x, Stress =="MeJA"))

SA_ps_l<-
  lapply(bac_alltaxa_phylosig_output_PCoA_nonsig$ordination_coordinates_l, function (x)
    subset_samples(x, Stress =="SA"))

# run the PoA axis pct and plot, bac
PCoA_528_plots_axis_pct_Control<-
  explained_PCoA_Axis_per_taxonomic_group(ps_obj_l = Control_ps_l)
names(PCoA_528_plots_axis_pct_Control)<-paste0(names(PCoA_528_plots_axis_pct_Control), "_Control")

PCoA_528_plots_axis_pct_MeJA<-
  explained_PCoA_Axis_per_taxonomic_group(ps_obj_l = MeJA_ps_l)
names(PCoA_528_plots_axis_pct_MeJA)<-paste0(names(PCoA_528_plots_axis_pct_Control), "_MeJA")

PCoA_528_plots_axis_pct_SA<-
  explained_PCoA_Axis_per_taxonomic_group(ps_obj_l = SA_ps_l)
names(PCoA_528_plots_axis_pct_SA)<-paste0(names(PCoA_528_plots_axis_pct_Control), "_SA")


# combine lists and save, bac
PCoA_528_3x_plots_axis_pct<-c(PCoA_528_plots_axis_pct_Control,
                              PCoA_528_plots_axis_pct_MeJA,
                              PCoA_528_plots_axis_pct_SA)

save(PCoA_528_3x_plots_axis_pct, file = "./14_Phylogenetic_distance/r_objects/PCoA_528_3x_plots_axis_pct.RData")






#split by stress, re-run function, Fun
Control_Fun_ps_l<-
  lapply(fun_alltaxa_phylosig_PCoA_nonsig$ordination_coordinates_l, function (x)
    subset_samples(x, Stress =="Control"))

MeJA_Fun_ps_l<-
  lapply(fun_alltaxa_phylosig_PCoA_nonsig$ordination_coordinates_l, function (x)
    subset_samples(x, Stress =="MeJA"))

SA_Fun_ps_l<-
  lapply(fun_alltaxa_phylosig_PCoA_nonsig$ordination_coordinates_l, function (x)
    subset_samples(x, Stress =="SA"))

# run the PoA axis pct and plot, Fun
PCoA_Fun_528_plots_axis_pct_Control<-
  explained_PCoA_Axis_per_taxonomic_group(ps_obj_l = Control_Fun_ps_l)
names(PCoA_Fun_528_plots_axis_pct_Control)<-paste0(names(PCoA_Fun_528_plots_axis_pct_Control), "_Control")

PCoA_Fun_528_plots_axis_pct_MeJA<-
  explained_PCoA_Axis_per_taxonomic_group(ps_obj_l = MeJA_Fun_ps_l)
names(PCoA_Fun_528_plots_axis_pct_MeJA)<-paste0(names(PCoA_Fun_528_plots_axis_pct_MeJA), "_MeJA")

PCoA_Fun_528_plots_axis_pct_SA<-
  explained_PCoA_Axis_per_taxonomic_group(ps_obj_l = SA_Fun_ps_l)
names(PCoA_Fun_528_plots_axis_pct_SA)<-paste0(names(PCoA_Fun_528_plots_axis_pct_SA), "_SA")

# combine lists and save, Fun
PCoA_Fun_528_3x_plots_axis_pct<-c(PCoA_Fun_528_plots_axis_pct_Control,
                              PCoA_Fun_528_plots_axis_pct_MeJA,
                              PCoA_Fun_528_plots_axis_pct_SA)

save(PCoA_Fun_528_3x_plots_axis_pct, file = "./14_Phylogenetic_distance/r_objects/PCoA_Fun_528_3x_plots_axis_pct.RData")



