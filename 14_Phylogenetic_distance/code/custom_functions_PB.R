# this script constains several custom functions to match distances in plant phylogeneis to distances in microbial communities for the rbassicacae family experiment
# most of these functions were written by Marcela Aragon
# some of these functions were written by Roland Berdaguer
# a few of these functions were sown together and streamlined by Pedro Beschoren


# /////////////////////////////////#
########## prepare_ps #############
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#

# prepare the ps object for the phylogeny pipeline (adjsut factors, remove sp outside the phylogenetic tree)
prepare_ps <- function(ps_obj) {

  # This function will adjust the factors in the metadata of a ps objet, and then...
  # it will remove, from the ps_objec, plant species that are not present in Bra_tree$tip.label (hard-coded and aivalable in global enviroment), and then..
  # i will remove ASVs with total zero reads after the filtering, and then...
  # it prints a few checks to see if the plant species retained are the ones you expect
  # the input is a phyloseq object
  # the output is a phyloseq object ready to continue wth the phylogeny pipeline

  ### NOTE: this function relies in metadata column positions to work. if those change, the function will crash
  ### NOTE: this function relies in a vector os plant_sp_abb present in the global enviroment as Bra_tree$tip.label. if it is missing, it will crash

 
  
  # Change some ps metadata to factor
  sample_data(ps_obj)[, c(1:9, 12:15, 17:25, 29:31)] <- lapply(sample_data(ps_obj)[, c(1:9, 12:15, 17:25, 29:31)], as.factor)

  # subistitute "_" for "." in Bo_F/M/S so it does not crash with next functions
  ps_obj@sam_data$Sp_abb_name <- gsub(pattern = "_", replacement = ".", x = ps_obj@sam_data$Sp_abb_name)

  # subset ps to only those spp found in the tree
  ps_sel <- subset_samples(ps_obj, Sp_abb_name %in% Bra_tree$tip.label)
  ps_sel <- prune_taxa(taxa_sums(ps_sel) > 0, ps_sel) # remove ASVs no longer present

  # just to check
  df_samples <- as.data.frame(with(sample_data(ps_sel), table(Sp_abb_name, Stress)))

  # print some check
  print("count of plant sp per stress in ps object")
  print(df_samples)
  print("are all plant sp in the bra_tree object? all shhould be TRUE!")
  print(ps_sel@sam_data$Sp_abb_name %>% unique() %in% Bra_tree$tip.label)
  print(Bra_tree$tip.label %in% ps_sel@sam_data$Sp_abb_name %>% unique())


  return(ps_sel)
}

########## DONE! #








# /////////////////////////////////#
## ps_to_distance_per_stress ######
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#


# define a function to calculate bray-curts and jaccard distances, from a list of BacFun ps objects, and returning the distance matrix split by stress
ps_to_distance_per_stress <- function(ps_BacFun_sel) {

  # this function takes a phyloseq object preapted by the prepare_ps() function
  # this input is a list of N phyloseq objects
  # it will calculate bray-curtis and jaccard dissimilarities for each of the 3 stresses in the brassicaceae family experiment
  # this function defines, inside itself, other functions that it is going to call
  # the output is a list of lists of 2 distances of  Bac or Fun communities, split by stress



  # Get distance matrix from Bacterial community, takes 5 min on highfilt data

  dist.bray <- lapply(ps_BacFun_sel, function(x) {
    phyloseq::distance(x, method = "bray")
  })
  dist.jaccard <- lapply(ps_BacFun_sel, function(x) {
    phyloseq::distance(x, method = "jaccard")
  })


  # make a list to use lapply afterwards, and also to save it as an RData object
  dist.list <- list(bray = dist.bray, jaccard = dist.jaccard)



  # Prepare matrix for mantel test using Roland's function (adapted)
  # which takes distance within species as well (diagonal is not 0)
  distance_to_Mantel <- function(dist, ps_l) {



    # This function takes a distance object and makes it a distance matrix ready

    # change object to matrix
    dist <- as.matrix(dist)

    # change object to data frame
    dist <- as.data.frame(dist)

    # row names to column "sample1"
    dist <- tibble::rownames_to_column(dist, "sample1")

    # transform wide to long
    dist <- dist %>% pivot_longer(
      cols = colnames(dist)[-1],
      names_to = "sample2",
      values_to = "distance"
    )

    ## replace sample names with plant species-stress treatment combination

    # get simplified metadata table with sample names, species codes and treatments
    metadata <- sample_data(ps_l)
    metadata <- as.matrix(metadata)
    metadata <- as.data.frame(metadata)
    # metadata$Sp_abb_name[metadata$Sp_abb_name ==  "Bo_M" ] <- 'Bo' #changing this so no problem with the "_" from Bo_M
    metadata <- metadata[, c("Sp_abb_name", "Stress")]
    metadata <- tibble::rownames_to_column(metadata, "sample_ID")
    metadata$species_treatment <- paste0(metadata$Sp_abb_name, "_", metadata$Stress)
    metadata <- metadata[, c("sample_ID", "species_treatment")]

    # merge dist with metadata
    dist <- merge(dist, metadata, by.x = "sample1", by.y = "sample_ID")
    dist$sample1 <- dist$species_treatment
    dist <- dist[, -4]
    dist <- merge(dist, metadata, by.x = "sample2", by.y = "sample_ID")
    dist$sample2 <- dist$species_treatment
    dist <- dist[, -4]

    # check how many rows have 0 in pairwise comparison (1,031 data points)
    nrow(dist[dist$distance == 0, ]) # check how many data points

    t <- as.data.frame(with(dist, table(sample1, sample2))) # 4,356 pairwise comparisons
    # lowest has 182 data points and highest has 256 (16*16)

    # calculate mean distances with the function "summarize"
    dist <- dist %>%
      filter(distance != 0) %>% # remove all those that have 0 (self comparison by sample)
      group_by(sample1, sample2) %>%
      dplyr::summarize(mean_distance = mean(distance))

    # [RB NOTE]: there is also a distance of a certain species/treatment combination with itself, which is the average distance between the replicates. we will later exclude the self-comparisons, but it is also good to check how the correlation looks if we include them.

    # subset by combination of sp and treatment
    dist <- separate(data = dist, col = sample1, into = c("species_code_1", "treatment_1"), sep = "_")
    dist <- separate(data = dist, col = sample2, into = c("species_code_2", "treatment_2"), sep = "_")

    # subset only rows that have Control for both
    dist_Control <- subset(dist, (treatment_1 == "Control" & treatment_2 == "Control"))
    # subset only rows that have MeJA for both
    dist_MeJA <- subset(dist, (treatment_1 == "MeJA" & treatment_2 == "MeJA"))
    # subset only rows that have SA for both
    dist_SA <- subset(dist, (treatment_1 == "SA" & treatment_2 == "SA"))

    # drop columns "treatment_1" and "treatment_2"
    dist_Control <- subset(dist_Control, select = -c(treatment_1, treatment_2))
    dist_MeJA <- subset(dist_MeJA, select = -c(treatment_1, treatment_2))
    dist_SA <- subset(dist_SA, select = -c(treatment_1, treatment_2))

    # Call function that uses "dist_X", which can then be either "dist_Control" or "dist_MeJA" or "dist_SA"
    # to make it the final distance matrix (to be used later)

    format_dist <- function(dist_X) {

      # transform long to wide [MA] important
      dist_X <- dist_X %>% pivot_wider(
        names_from = "species_code_2",
        values_from = "mean_distance"
      )

      # change to matrix, assign row names
      dist_X <- as.matrix(dist_X)
      rownames(dist_X) <- dist_X[, 1] # change row 1 to row names
      dist_X <- dist_X[, -1] # remove the column that contained the row names

      return(dist_X) # [MA] returns distance matrix
    }


    # run the function on the Stress- bray-curtis distance matrices
    dist_Control <- format_dist(dist_Control)
    dist_MeJA <- format_dist(dist_MeJA)
    dist_SA <- format_dist(dist_SA)


    # make it into a list and return it out of the function
    dist_ByStress_list <- list(dist_Control, dist_MeJA, dist_SA)
    names(dist_ByStress_list) <- c("Control", "MeJA", "SA")

    # change to numeric
    numeric_dist_ByStress_list <- lapply(dist_ByStress_list, function(x) {
      class(x) <- "numeric"
      return(x)
    })

    return(numeric_dist_ByStress_list)
  }



  # apply function
  dist.bray_l <- mapply(function(x, y) {
    distance_to_Mantel(dist = x, ps_l = y)
  },
  x = dist.list$bray,
  y = ps_BacFun_sel,
  SIMPLIFY = FALSE
  )



  dist.jaccard_l <- mapply(function(x, y) {
    distance_to_Mantel(dist = x, ps_l = y)
  },
  x = dist.list$jaccard,
  y = ps_BacFun_sel,
  SIMPLIFY = FALSE
  )


  output <- list(
    "dist.bray_l" = dist.bray_l,
    "dist.jaccard_l" = dist.jaccard_l
  )

  return(output)
}

########## DONE! #





























# /////////////////////////////////#
### matrix_order_By_phylo #########
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#

# function to change the order in matrix
matrix_order_By_phylo <- function(x) {
  # this function will change the order of the species in the distnce matrix to amthc the phylogenetic tree
  # it's input is the output of the ps_to_distance_per_stress() function
  # it's output is a reordered version of the ps_to_distance_per_stress() output
  # NOTE: it relies on the object "phylo_tips", which is defined in the global environment
  phylo_tips <- Bra_tree$tip.label
  
  new_order_matrix <- x[phylo_tips, phylo_tips]
  return(new_order_matrix)
}


########## DONE! #



# /////////////////////////////////#
########## check_order#############
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#

# function to check order of species is okay
check_order <- function(x) {
  # this function will check if the order of plant species in the disntances of plant phylogenies matches the distances of the microbial communities
  # it's input is the output of the matrix_order_By_phylo() function
  # It's output is a series of logical values, which should all be true if the species order amtch. they should all be mathicng to continue!
  # NOTE: this fucntion relies on the object pd, which are the distances of the plant species in the base substitution trees
  print("Are the names in the plant distances in the same order as the names in the microbial community distances?")
  print(rownames(pd) == rownames(x))
}

########## DONE! #























# /////////////////////////////////#
########## df_Mantel ##############
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#


# extract stats from a mantel test
df_Mantel <- function(x) {

  # this function will extract p, R, and treatment emtadata values from a matel test, and then
  # it saves them as a data frame
  # the function is called inside distance_per_stress_to_mantel_df()
  # the input (mantel tests per plant stress) is provided by distance_per_stress_to_mantel_df()

  Mantel_R <- as.data.frame(x$statistic)
  pvalue <- as.data.frame(x$signif)
  Treatment <- deparse(substitute(x)) # gets name of the object as a string character

  df <- data.frame(Treatment, Mantel_R, pvalue)

  return(df)
}

########## DONE! #


# /////////////////////////////////#
# distance_per_stress_to_mantel_df ####
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#

# define a function to calculate mantel tests and then extract p, r values from it
distance_per_stress_to_mantel_df <- function(distance_per_stress_output, correlation_method, plant_distances = pd) {

  # this function will calculate mantel distances between plant phylogenies and microbial communities (both bray-curtis and jaccard), and then...
  # it will save the key metric output into a df
  # the input for this function is the output of the ps_to_distance_per_stress() function
  # the output of this function is a df with matel test results
  # you can also define the correlation tests for mante, and add a specific plant_distances object

  # Calculate Mantel test and make it a df for each distance matrix

  # Bray-Curtis
  mantel.bray <- lapply_lapply(distance_per_stress_output$dist.bray_l, function(x) {
    set.seed(123)
    vegan::mantel(x, plant_distances, method = correlation_method, na.rm = FALSE, permutations = 999)
  })

  # make df for each element
  df.Mantel_Bray <- lapply_lapply(mantel.bray, df_Mantel)
  
  print("CAUTION: if the name of your ps_object_l includes a dot, this function will fail!")

  # merging df
  df.Mantel_Bray <- as.data.frame(unlist(df.Mantel_Bray)) %>%
    tibble::rownames_to_column() %>%
    tidyr::separate(
      col = rowname, into = c("Community", "Treatment", "metric"),
      extra = "merge",
      sep ="\\."
    ) %>% # takes everything after first point and puts it together
    mutate(distance = "Bray") %>% # change this manually for distance
    subset(metric != "Treatment") %>%
    dplyr::rename_with(.cols = 4, ~"value")

  # Jaccard
  mantel.jaccard <- lapply_lapply(distance_per_stress_output$dist.jaccard_l, function(x) {
    set.seed(123)
    vegan::mantel(x, plant_distances, method = correlation_method, na.rm = FALSE, permutations = 999)
  })

  # make df for each element
  df.Mantel_Jaccard <- lapply_lapply(mantel.jaccard, df_Mantel)

  # merging df
  df.Mantel_Jaccard <- as.data.frame(unlist(df.Mantel_Jaccard)) %>%
    tibble::rownames_to_column() %>%
    tidyr::separate(
      col = rowname, into = c("Community", "Treatment", "metric"),
      extra = "merge",
      sep ="\\."
    ) %>% # takes everything after first point and puts it together
    mutate(distance = "Jaccard") %>% # change this manually for distance
    subset(metric != "Treatment") %>%
    dplyr::rename_with(.cols = 4, ~"value")


  # Merge df
  df.Mantel <- rbind(df.Mantel_Bray, df.Mantel_Jaccard)

  # make it wide
  df.Mantel <- df.Mantel %>%
    distinct() %>% # removes duplicate values
    pivot_wider(names_from = metric, values_from = value)

  # change colnames
  df.Mantel <- df.Mantel %>%
    dplyr::rename_with(.cols = 4, ~"r2_Mantel") %>%
    dplyr::rename_with(.cols = 5, ~"p_value") %>%
    dplyr::mutate_at(c("r2_Mantel", "p_value"), as.numeric) %>%
    dplyr::mutate_at(c("Treatment", "distance"), as.factor) %>%
    dplyr::mutate(significant = dplyr::case_when(
      p_value < 0.05 ~ "Yes",
      p_value > 0.05 ~ "No"
    ))
  # round up R2
  df.Mantel$r2_Mantel <- round(df.Mantel$r2_Mantel, digit = 3)

  return(df.Mantel)
}

########## DONE! #




# /////////////////////////////////#
# a few ggplot paramethers... ####
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#

# Set axis looks
theme_set(theme_bw())
axis_looks <- theme(axis.text.x = element_text(
  colour = "black", size = 10,
  face = "bold"
)) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 15, face = "bold"))

library(RColorBrewer)
ColorDistance <- brewer.pal(4, "BrBG") # to get 4 contrasting colors

########## DONE! #











# /////////////////////////////////#
#### plot_mantel_correlations #####
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#

# this function will create a plot showing mantel test results
plot_mantel_correlations <- function(distance_per_stress_to_mantel_df_output, quoted_plot_title) {

  # this function takes the output of the distance_per_stress_to_mantel_df() function to make plots of p and r values suitable to the brassicaceae family experiment design
  # the input is the utput of distance_per_stress_to_mantel_df()
  # the output is a pannel of ggplots

  # plot p values
  p.summary <- ggplot(
    data = distance_per_stress_to_mantel_df_output,
    mapping = aes(x = Treatment, y = p_value, fill = distance)
  ) +
    geom_point(aes(colour = distance),
      position = position_dodge(0.8), size = 3
    ) +
    scale_colour_manual(values = ColorDistance)

  plot_summary <- p.summary + 
    facet_wrap(~Community) +
    geom_hline(aes(yintercept = 0.05), colour = "red", size = 1) +
    axis_looks +
    theme(legend.position = "top")


  # plot R2 values
  r.summary <- ggplot(
    data = distance_per_stress_to_mantel_df_output,
    mapping = aes(x = Treatment, y = r2_Mantel, fill = distance)
  ) +
    geom_point(aes(colour = distance),
      position = position_dodge(0.8), size = 3
    ) +
    scale_colour_manual(values = ColorDistance) +
    axis_looks

  plot_r.summary <- r.summary + facet_wrap(~Community) +
    geom_hline(aes(yintercept = 0.0), colour = "black", size = 1) +
    axis_looks +
    theme(legend.position = "top")

  output <- ggarrange(plot_summary, plot_r.summary, labels = quoted_plot_title)

  return(output)
}

########## DONE! #






# define function to turn mantel correlations (output of distance_per_stress()) into a df for aplot
mantel_correlation_raw<- function(dist_X){
  # dist_X = output of distance_per_stress function
  #this function was done by Roland Berdaguer  
  
  # pivot wide to long
  pd_df <- as.data.frame(pd) #pd are the plant distances, defined n the local enviroment from the Bra_tree
  pd_df$species1 <- rownames(pd_df)
  pd_df_long <- pd_df %>% pivot_longer(cols=colnames(pd_df)[-ncol(pd_df)], # we exclude the last column  
                                       names_to='species2',                  # because contains'species1'
                                       values_to='genetic_distance')
  
  
  dist_X_df <- as.data.frame(dist_X)
  dist_X_df$species1 <- rownames(dist_X_df)
  dist_X_df_long <- dist_X_df %>% pivot_longer(cols=colnames(dist_X_df)[-ncol(dist_X_df)], # we exclude the last column because it contains "species1"
                                               names_to='species2',
                                               values_to='microbiome_distance')
  
  # merge tables
  merged_distances <- merge(pd_df_long, dist_X_df_long, by = c("species1", "species2"))
  
  # remove inverted duplicate rows (e.g. At-Slm and Slm-At)
  merged_distances <- merged_distances %>%
    group_by(grp = paste(pmax(species1, species2), pmin(species1, species2), sep = "_")) %>%
    dplyr::slice(1) #[MA] what is this doing? 
  
  #pmax puts always the first sp by alphabetical order, so it will give the same order under the 
  #'grp' column, for instance Al-Bo and Bo-Al will be both Al_Bo.
  #With 'slice(1)' then you are keeping the first entry of every unique group on the 'grp'
  #column and hence removing the repeated ones. 
  
  
  return(merged_distances)
  
}

# define function to assign pairwise phylogeny comparisons
pairwise_comparison_byLineage <- function(x){
  
  lineageI <- c("At", "Al", "Cas", "Tg", "Ec", "Mm", "Bv",
                "Lr", "Ds")
  
  lineageII <- c("Bi", "Bn", "Sar", "Hi", "Br", "Bo.M", "Dt", "Si", "It") 
  ext_lineageII<-  c("Co", "Lm", "Ia") 
  
  lineageIII <-  c("Es")
  
  #make df with a new column in which categories are allocated 
  df  <-     x %>% 
    dplyr::mutate(sp1_lineage=dplyr::case_when( #classification of sp1
      species1 %in% lineageII  ~ "lineage_2",
      species1 %in% lineageI  ~ "lineage_1",
      species1 %in% lineageIII  ~ "lineage_3",
      species1 %in% ext_lineageII  ~ "ext_lineage_2",
      TRUE ~ "Other")) %>% 
    dplyr::mutate(sp2_lineage=dplyr::case_when( #classification of sp2
      species2 %in% lineageII  ~ "lineage_2",
      species2 %in% lineageI  ~ "lineage_1",
      species2 %in% lineageIII  ~ "lineage_3",
      species2 %in% ext_lineageII  ~ "ext_lineage_2",
      TRUE ~ "Other")) %>% 
    dplyr::mutate(pairwise_comparison=dplyr::case_when( #pairwise_contrasts
      sp1_lineage == "lineage_1" & sp2_lineage == "lineage_1" ~ "self-I",
      sp1_lineage == "lineage_2" & sp2_lineage == "lineage_2" ~ "self-II",
      sp1_lineage == "lineage_3" & sp2_lineage == "lineage_3" ~ "self-III",
      sp1_lineage == "ext_lineage_2" & sp2_lineage == "ext_lineage_2" ~ "self-ext_II",
      
      ((sp1_lineage == "lineage_1" & sp2_lineage == "lineage_2") |
         (sp1_lineage == "lineage_2" & sp2_lineage == "lineage_1")) ~ "I-II",
      ((sp1_lineage == "lineage_1" & sp2_lineage == "lineage_3") |
         (sp1_lineage == "lineage_3" & sp2_lineage == "lineage_1")) ~ "I-III",
      ((sp1_lineage == "lineage_2" & sp2_lineage == "lineage_3") |
         (sp1_lineage == "lineage_3" & sp2_lineage == "lineage_2")) ~ "II-III",
      ((sp1_lineage == "lineage_1" & sp2_lineage == "ext_lineage_2") |
         (sp1_lineage == "ext_lineage_2" & sp2_lineage == "lineage_1")) ~ "I-ext_II",
      ((sp1_lineage == "lineage_2" & sp2_lineage == "ext_lineage_2") |
         (sp1_lineage == "ext_lineage_2" & sp2_lineage == "lineage_2")) ~ "II-ext_II",
      ((sp1_lineage == "lineage_3" & sp2_lineage == "ext_lineage_2") |
         (sp1_lineage == "ext_lineage_2" & sp2_lineage == "lineage_3")) ~ "III-ext_II",
      TRUE ~ "Other")) #assigns 'Other' to anything else that's not on those categories (just to check)
  
  #change new columns as factor
  df[c(1:2,8)] <- lapply(df[c(1:2,8)], factor)
  
  #order by distance lower to higher
  df$pairwise_comparison <- factor(df$pairwise_comparison, 
                                   levels = c("self-I", 
                                              "self-II", 
                                              "self-III", 
                                              "self-ext_II", 
                                              "II-III", 
                                              "I-III", 
                                              "I-II",
                                              "I-ext_II",
                                              "II-ext_II",
                                              "III-ext_II"))
  
  print(df %>% 
          dplyr::group_by(pairwise_comparison) %>% 
          dplyr::count())#check
  
  return(df)
}

# define function to plot pairwise phylogeny comparisons
plot_pairwise_lineage <-  function (x){
  # x = one list element of the pairwise_comparison_byLineage() output
  
  #make vector with colors (Viridis, inferno + ranom colors)
  color_pairwise_lineage <- c("#000004","#420a68", "#932667", "#dd513a", "#fca50a", "#fcffa4","#D55E00","#0072B2","#CC79A7","#56B4E9")
  
  #some plot paramenters
  titlenames <- c("Ctr", "MeJA", "SA")
  
  three_plots<-purrr::map2(x, titlenames, function(x, titlenames) 
    ggplot(x, aes(x=genetic_distance, y=microbiome_distance)) +
      geom_point(size=7,alpha=0.7, shape=21, stroke=0.9,
                 aes(fill=pairwise_comparison), color='black')+
      scale_fill_manual(values=color_pairwise_lineage)+
      scale_color_manual(values=color_pairwise_lineage)+
      xlab("Plant phylogenetic distance")+
      ylab("community Bray-Curtis dissimilarity")+
      geom_smooth(method='lm', color="black")+
      ggtitle(titlenames)+
      axis_looks)
  
  output<-ggarrange(three_plots$Control,
                    three_plots$MeJA,
                    three_plots$SA,
                    ncol = 1,
                    labels = "AUTO",
                    widths = 20,
                    heights = 20,
                    common.legend = TRUE,
                    legend = "top",
                    unit = "cm")
  return(output)
}




# /////////////////////////////////#
######## ps_to_mantel_pipeline ########
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#


# define a function to run the entire pipeline of custom fucntions above
ps_to_mantel_pipeline <- function(ps_obj_l, ...) {
  # this function will run the  phylogenetic pipeline that goes from ps object to a matel plot
  # the input is a phyloseq object
  # the output is a mantel test of the community in the ps object agasint the plant phylogenies in the tree

  
  
  # laod and prepare Bra_tree. this fucntion anme Brassica oleraceae as Bo_M
  Bra_tree <- load_prepare_MY_tree()
  Bra_tree$tip.label[ Bra_tree$tip.label == "Bo_M"]<-"Bo.M"
  phylo_tips <- Bra_tree$tip.label
  # overwrite phylo_tips becuase the next fucntions reli on Bo.M and not on Bo_M
#  phylo_tips <- c(
 #   "At", "Al", "Cas", "Tg", "Ec", "Mm", "Bv",
  #  "Lr", "Ds", "Bi", "Bn", "Sar", "Hi", "Br", "Bo.M", # changed Bo_M to Bo so it matches code below
   # "Dt", "Si", "It", "Co", "Lm", "Ia", "Es"
  #)

  # change for new names
  #Bra_tree$tip.label <- phylo_tips
  
  #check bra tree
#  print(Bra_tree)
  print(Bra_tree$tip.label)
  

  # calculate plant phylogenetic distance Fitzpatrick et al., 2019
  pd <- cophenetic.phylo(Bra_tree)
  
  #save "Bra_tree" and "pd" in the global enviroment so we can make some plots
  assign("Bra_tree", Bra_tree, envir=globalenv())
  assign("pd", pd, envir=globalenv())


  # prepare phyloseq object for the rest of the phylogeny pipeline (subset samples, adjust factors)
  ps_obj_l_ready <- lapply(ps_obj_l, function(x) prepare_ps(ps_obj = x))

  # define bray-curtis and jaccard distances for Bac anf Fun under Control, MeJA and SA conditions with a custom function
  distance_per_stress <- ps_to_distance_per_stress(ps_obj_l_ready)

  # change the order os the plant species in the bray-curties tree to reflect the phylogeny tree
  distance_per_stress <- lapply(distance_per_stress, function(x) lapply_lapply(x, matrix_order_By_phylo))

  # check if orders match; they should all be TRUE!
  lapply(distance_per_stress, function(x) lapply_lapply(x, check_order))

  # run mantel tests of microbial community with plant phylogenies, and then return a df summarizing the results
  mantel_df_output <- distance_per_stress_to_mantel_df(distance_per_stress, correlation_method = "kendall", plant_distances = pd)

  # plot mantel correlations as in the custom function
  plot_output <- plot_mantel_correlations(mantel_df_output, quoted_plot_title = "")
  
  #apply function to get list of 3 df's to make plot afterwards 
  corr.Bray <- lapply_lapply(distance_per_stress$dist.bray_l, mantel_correlation_raw)
  
  #make a df with pairwise lineage comparsion data
  corr.Bray_df<-lapply_lapply(corr.Bray, pairwise_comparison_byLineage)
  
  # make a plot form this df
  bray_correlation_plot<-lapply(corr.Bray_df,plot_pairwise_lineage)
  
  

  output<-list("mantel_df_output" = mantel_df_output,
               "plot_output" = plot_output,
               "bray_correlation_plot" = bray_correlation_plot)
  


  # ends
  return(output)
}

########## DONE! #



# /////////////////////////////////#
############ ps_to_df #############
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#


# turns a single ps bject itno a df ready for phylogenetic signal analysis
ps_to_df <- function(ps_obj, colstokeep) {
  # this function will take a ps object from the rbassicaceae family experiment and turn it into a subset metadata dataframe separated by stress, in a list
  # it has two inputs: one is a a phyloseq object
  # the second input is a vector of metadata column names to keep and test against the phylogeny
  # the output is a dataframe of sample emtadata suitable to be tested agaisnt the phylogeny


  df <- as.data.frame(sample_data(ps_obj))
  colnames(df)
  colstokeep <- colstokeep # c("Sp_abb_name", "Stress", "root_weight", "Shoot_.DryWeight", "Host_DNA_contamination_pct", "library_size") # pedro removed a few metadata columns that would not help ;  more per-species numeric data can be added here (n of DA ASVs, permaova R, PC1 position...)

  df_plant <- df[, (colnames(df) %in% colstokeep)]
  colnames(df_plant)
  df_plant <- as_tibble(df_plant) # needs to be change as.tibble otherwise summarise doesn't work
  str(df_plant)

  # change name so it doesn't mess up with the underscore
  # levels(df_plant$Sp_abb_name)[match("Bo_M",levels(df_plant$Sp_abb_name))] <- "Bo"

  # Silenced in 24/March/23: we don't need root_weight as an obligaory variable on every tree (inflates global errors)
  #  df_plant <- as.data.frame(df_plant) #change back to dataframe to avoid errors
  #  df_plant$root_weight <- as.numeric(df_plant$root_weight ) #change to numeric
  #  str(df_plant) #check

  # make new df with mean to be used for phylogenetic signal
  # [MA] I still need to get also the 'N' value to calculate standard error but I couldn't do it!!
  df_plant <- df_plant %>%
    dplyr::group_by(Sp_abb_name, Stress) %>%
    dplyr::summarise(dplyr::across(everything(), list(mean = mean), na.rm = TRUE)) %>%
    dplyr::relocate(Sp_abb_name)

  # order by phylogeny
  df_plant$Sp_abb_name <- factor(df_plant$Sp_abb_name,
    levels = levels(ps_obj@sam_data$Sp_abb_name)
  ) # changed by pedro, was targeting unexisting object

  # make 3 df, one for each stress
  df_plant_Control <- df_plant[df_plant$Stress == "Control", ]
  df_plant_MeJA <- df_plant[df_plant$Stress == "MeJA", ]
  df_plant_SA <- df_plant[df_plant$Stress == "SA", ]

  # make it a list to use lapply after
  df_plant_l <- list(df_plant_Control, df_plant_MeJA, df_plant_SA)
  names(df_plant_l) <- c("Control", "MeJA", "SA")

  return(df_plant_l)
}

########## DONE! #













# ////////////////////////////////////////#
#### plantdf_to_Abouheif_and_moran ####
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#

# define a function to take the prepared df and calculate  Abouheif_and_moran agasint a phylogenetic tree on it
plantdf_to_Abouheif_and_moran <- function(ps_to_df_output, phylogenetic_tree, show_abouheif_histogram = FALSE) {
  # This funtion will take a df in the format of the ps_to_df() fucntion and match such emtadata agasint the phylogenetic tree
  # it has two arguments: the output of ps_to_df() and a phylogenetic tree like Bra_tree
  # the output are the metrics for the phylogenetic signal (Abouheif and abouheif )


  # prepare data
  metadata_for_phylogeny <- as.data.frame(ps_to_df_output, row.names = NULL) # change to df, otherwise weird errors

  row.names(metadata_for_phylogeny) <- metadata_for_phylogeny[, 1] # add rownames based on sp column
  metadata_for_phylogeny <- metadata_for_phylogeny[, -(1:2)] # remove th two first columns 'Sp' and "Stress" column


  # Moran's
  phylotraits <- phylo4d(phylogenetic_tree, metadata_for_phylogeny) # from phylobase package
  # Then, we do moran test using some Monte Carlo simulations (default is 999 randomizations)
  set.seed(12345)
  moran.test <- abouheif.moran(phylotraits, method = "Abouheif", nrepet = 9999)

  # Albouheif's
  set.seed(12345)
  abouheif.test <- abouheif.moran(phylotraits, method = "oriAbouheif", nrepet = 9999)



  output <- list(
    "moran.test" = moran.test,
    "abouheif.test" = abouheif.test
  )

  #plots showing the measured phylogenetic signal agasint the 999 permutations is silenced for now are being sillenced
 # plot(moran.test, main = "")
#  title(main = "moran.test")
  if (show_abouheif_histogram == TRUE) { # this is an arugment of phylogenetic_signal_pipeline()

  plot(abouheif.test, main = "")
  title(main = "abouheif.test")
  
  }

  return(output)
}

########## DONE! #









# /////////////////////////////////#
##### phylogenetic_lambda_K #######
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#


# dcalculate phylogenetic signal on a trait present in the previous df
phylogenetic_lambda_K <- function(ps_to_df_output, selected_trait, phylogenetic_tree) {
  # this function will calculate lambda and K of a trait against the phylogeny
  # it has 3 inputs: the output of the ps_to_df() function, the quoted name of metadata column present in that df, and a phylogenetic tree like Bra_tree
  # the output si a list with K and lambda values for the trait

  # prepare data
  test <- as.data.frame(ps_to_df_output, row.names = NULL) # change to df, otherwise weird errors
  # test <- test[,-c(1)]
  test <- test %>%
    dplyr::relocate(Sp_abb_name)

  row.names(test) <- test[, 1] # add rownames
  test <- test[, -1] # remove 'Sp' column
  # test <- test[,grepl("mean",names(test))] #keep only mean values for now




  # test data
  trait <- test[, selected_trait] # shoot dry weight
  names(trait) <- rownames(test)

  # We choose plant biomass as the trait we are testing for phylogenetic signal. Then, we do the test with 999 randomizations:
  set.seed(12345)
  lambda <- phylosig(phylogenetic_tree, trait, method = "lambda", test = TRUE, nsim = 9999) # sig!

  # From Munkermuller et al., 2012: stronger deviations from zero indicate stronger rela- tionships between trait values and the phylogeny, a value close to zero indicates phylogenetic independence and a value of one indicates that species’ traits are distributed as expected under BM.In most cases, the upper limit of Pagel’s k is close to one (see Materials and methods for details), while Blomberg’s K can take higher values indicating stronger trait similarity between related speciesthanexpectedunder BM

  # Bloomberg's
  # same settings as for Page's l
  set.seed(12345)
  k <- phylosig(phylogenetic_tree, trait, method = "K", test = TRUE, nsim = 9999) # sig!

  output <- list(
    "lambda" = lambda,
    "k" = k
  )

  return(output)
}





#function t extract results from Pagel's lambda

lambda_k_results <- function(input_l){
  
  lambda_value_l <- lapply(input_l, function(x) x$lambda$lambda)
  lambda_pvalue_l <- lapply(input_l, function(x) x$lambda$P)
  k_value_l <- lapply(input_l, function(x) x$k$K)
  k_pvalue_l <- lapply(input_l, function(x) x$k$P)
  
  # Convert the list to a named vector
  lambda_value_v <- unlist(lambda_value_l, recursive = FALSE)
  lambda_pvalue_v <- unlist(lambda_pvalue_l, recursive = FALSE)
  k_value_v <- unlist(k_value_l, recursive = FALSE)
  k_pvalue_v <- unlist(k_pvalue_l, recursive = FALSE)
  
  # Stack the vector and create a data frame
  lambda_value_df <- as.data.frame(stack(lambda_value_v))
  lambda_pvalue_df <- as.data.frame(stack(lambda_pvalue_v))
  lambda_df <- merge(lambda_value_df, lambda_pvalue_df, by = "ind")
  
  colnames(lambda_df) <- c("stress", "value", "pvalue") 
  
  lambda_df <- mutate(lambda_df, index = "Pagel's lambda") #specify index
  
  k_value_df <- as.data.frame(stack(k_value_v))
  k_pvalue_df <- as.data.frame(stack(k_pvalue_v))
  k_df <- merge(k_value_df, k_pvalue_df, by = "ind")
  
  colnames(k_df) <- c("stress", "value", "pvalue") 
  
  k_df <- mutate(k_df, index = "Blomberg's K")#specify index
  
  #merge both dataframes
  output <- rbind(lambda_df, k_df)}


#done!















####################### separating abcteria and fungi on the cutom functions


# /////////////////////////////////#
## ps_to_distance_per_stress_1kingdom ######
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#

# define a function to calculate bray-curts and jaccard distances, from a list of BacFun ps objects, and returning the distance matrix split by stress
ps_to_distance_per_stress_1kingdom <- function(ps_Bac_OR_Fun_sel) {

  # this function takes a phyloseq object preapted by the prepare_ps() function
  # NOTE: this input is a list of 2 phyloseq objects: one "Bac" and one "Fun". objects with other names will crash!
  # it will calculate bray-curtis and jaccard dissimilarities for each of the 3 stresses in the brassicaceae family experiment
  # this function defines, inside itself, other functions that it is going to call
  # the output is a list of lists of 2 distances of  Bac or Fun communities, split by stress



  # Get distance matrix from Bacterial community, takes 5 min on highfilt data

  dist.bray <-
    phyloseq::distance(ps_Bac_OR_Fun_sel, method = "bray")
  dist.jaccard <-
    phyloseq::distance(ps_Bac_OR_Fun_sel, method = "jaccard")


  # make a list to use lapply afterwards, and also to save it as an RData object
  dist.list <- list(bray = dist.bray, jaccard = dist.jaccard)



  # Prepare matrix for mantel test using Roland's function (adapted)
  # which takes distance within species as well (diagonal is not 0)
  distance_to_Mantel_Bac <- function(dist) {


    # dist <- dist.bray$Control

    # This function takes a distance object and makes it a distance matrix ready

    # change object to matrix
    dist <- as.matrix(dist)

    # change object to data frame
    dist <- as.data.frame(dist)

    # row names to column "sample1"
    dist <- tibble::rownames_to_column(dist, "sample1")

    # transform wide to long
    dist <- dist %>% pivot_longer(
      cols = colnames(dist)[-1],
      names_to = "sample2",
      values_to = "distance"
    )

    ## replace sample names with plant species-stress treatment combination

    # get simplified metadata table with sample names, species codes and treatments
    metadata <- sample_data(ps_Bac_OR_Fun_sel)
    metadata <- as.matrix(metadata)
    metadata <- as.data.frame(metadata)
    # metadata$Sp_abb_name[metadata$Sp_abb_name == "Bo_M"] <- 'Bo' #changing this so no problem with the "_" from Bo_M
    metadata <- metadata[, c("Sp_abb_name", "Stress")]
    metadata <- tibble::rownames_to_column(metadata, "sample_ID")
    metadata$species_treatment <- paste0(metadata$Sp_abb_name, "_", metadata$Stress)
    metadata <- metadata[, c("sample_ID", "species_treatment")]

    # merge dist with metadata
    dist <- merge(dist, metadata, by.x = "sample1", by.y = "sample_ID")
    dist$sample1 <- dist$species_treatment
    dist <- dist[, -4]
    dist <- merge(dist, metadata, by.x = "sample2", by.y = "sample_ID")
    dist$sample2 <- dist$species_treatment
    dist <- dist[, -4]

    # check how many rows have 0 in pairwise comparison (1,031 data points)
    nrow(dist[dist$distance == 0, ]) # check how many data points

    t <- as.data.frame(with(dist, table(sample1, sample2))) # 4,356 pairwise comparisons
    # lowest has 182 data points and highest has 256 (16*16)

    # calculate mean distances with the function "summarize"
    dist <- dist %>%
      filter(distance != 0) %>% # remove all those that have 0 (self comparison by sample)
      group_by(sample1, sample2) %>%
      dplyr::summarize(mean_distance = mean(distance))

    # [RB NOTE]: there is also a distance of a certain species/treatment combination with itself, which is the average distance between the replicates. we will later exclude the self-comparisons, but it is also good to check how the correlation looks if we include them.

    # subset by combination of sp and treatment
    dist <- separate(data = dist, col = sample1, into = c("species_code_1", "treatment_1"), sep = "_")
    dist <- separate(data = dist, col = sample2, into = c("species_code_2", "treatment_2"), sep = "_")

    # subset only rows that have Control for both
    dist_Control <- subset(dist, (treatment_1 == "Control" & treatment_2 == "Control"))
    # subset only rows that have MeJA for both
    dist_MeJA <- subset(dist, (treatment_1 == "MeJA" & treatment_2 == "MeJA"))
    # subset only rows that have SA for both
    dist_SA <- subset(dist, (treatment_1 == "SA" & treatment_2 == "SA"))

    # drop columns "treatment_1" and "treatment_2"
    dist_Control <- subset(dist_Control, select = -c(treatment_1, treatment_2))
    dist_MeJA <- subset(dist_MeJA, select = -c(treatment_1, treatment_2))
    dist_SA <- subset(dist_SA, select = -c(treatment_1, treatment_2))

    # Call function that uses "dist_X", which can then be either "dist_Control" or "dist_MeJA" or "dist_SA"
    # to make it the final distance matrix (to be used later)

    format_dist <- function(dist_X) {

      # transform long to wide [MA] important
      dist_X <- dist_X %>% pivot_wider(
        names_from = "species_code_2",
        values_from = "mean_distance"
      )

      # change to matrix, assign row names
      dist_X <- as.matrix(dist_X)
      rownames(dist_X) <- dist_X[, 1] # change row 1 to row names
      dist_X <- dist_X[, -1] # remove the column that contained the row names

      return(dist_X) # [MA] returns distance matrix
    }


    # run the function on the Stress- bray-curtis distance matrices
    dist_Control <- format_dist(dist_Control)
    dist_MeJA <- format_dist(dist_MeJA)
    dist_SA <- format_dist(dist_SA)


    # make it into a list and return it out of the function
    dist_ByStress_list <- list(dist_Control, dist_MeJA, dist_SA)
    names(dist_ByStress_list) <- c("Control", "MeJA", "SA")

    # change to numeric
    numeric_dist_ByStress_list <- lapply(dist_ByStress_list, function(x) {
      class(x) <- "numeric"
      return(x)
    })

    return(numeric_dist_ByStress_list)
  }



  # apply function
  dist.bray_l <- distance_to_Mantel_Bac(dist.list$bray)
  dist.jaccard_l <- distance_to_Mantel_Bac(dist.list$jaccard)

  output <- list(
    "dist.bray_l" = dist.bray_l,
    "dist.jaccard_l" = dist.jaccard_l
  )

  return(output)
}

########## DONE! #








# /////////////////////////////////#
# distance_per_stress_to_mantel_df_1kingdom ####
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#

# define a function to calculate mantel tests and then extract p, r values from it
distance_per_stress_to_mantel_df_1kingdom <- function(distance_per_stress_1kingdom_output, correlation_method) {

  # this function will calculate mantel distances between plant phylogenies and microbial communities (both bray-curtis and jaccard), and then...
  # it will save the key metric output into a df
  # the input for this function is the output of the ps_to_distance_per_stress() function
  # the output of this function is a df with matel test results

  # Calculate Mantel test and make it a df for each distance matrix

  # Bray-Curtis
  mantel.bray <- lapply(distance_per_stress_1kingdom_output$dist.bray_l, function(x) {
    set.seed(123)
    vegan::mantel(x, pd, method = correlation_method, na.rm = FALSE, permutations = 999)
  })

  # make df for each element
  df.Mantel_Bray <- lapply(mantel.bray, df_Mantel)

  # merging df
  df.Mantel_Bray <- as.data.frame(unlist(df.Mantel_Bray)) %>%
    tibble::rownames_to_column() %>%
    tidyr::separate(
      col = rowname, into = c("Treatment", "Community", "metric"),
      extra = "merge"
    ) %>% # takes everything after first point and puts it together
    mutate(distance = "Bray") %>% # change this manually for distance
    subset(metric != "Treatment") %>%
    dplyr::rename_with(.cols = 4, ~"value")

  df.Mantel_Bray <- df.Mantel_Bray[, -2]

  # Jaccard
  mantel.jaccard <- lapply(distance_per_stress_1kingdom_output$dist.jaccard_l, function(x) {
    set.seed(123)
    vegan::mantel(x, pd, method = correlation_method, na.rm = FALSE, permutations = 999)
  }) # Ctr=NS, MeJA=NS, SA=Sig

  # make df for each element
  df.Mantel_Jaccard <- lapply(mantel.jaccard, df_Mantel)

  # merging df
  df.Mantel_Jaccard <- as.data.frame(unlist(df.Mantel_Jaccard)) %>%
    tibble::rownames_to_column() %>%
    tidyr::separate(
      col = rowname, into = c("Treatment", "Community", "metric"),
      extra = "merge"
    ) %>% # takes everything after first point and puts it together
    mutate(distance = "Jaccard") %>% # change this manually for distance
    subset(metric != "Treatment") %>%
    dplyr::rename_with(.cols = 4, ~"value")

  df.Mantel_Jaccard <- df.Mantel_Jaccard[, -2]

  # Merge df
  df.Mantel <- rbind(df.Mantel_Bray, df.Mantel_Jaccard)

  # make it wide
  df.Mantel <- df.Mantel %>%
    distinct() %>% # removes duplicate values
    pivot_wider(names_from = metric, values_from = value)

  # change colnames
  df.Mantel <- df.Mantel %>%
    dplyr::rename_with(.cols = 3, ~"r2_Mantel") %>%
    dplyr::rename_with(.cols = 4, ~"p_value") %>%
    dplyr::mutate_at(c("r2_Mantel", "p_value"), as.numeric) %>%
    dplyr::mutate_at(c("Treatment", "distance"), as.factor) %>%
    dplyr::mutate(significant = dplyr::case_when(
      p_value < 0.05 ~ "Yes",
      p_value > 0.05 ~ "No"
    ))
  # round up R2
  df.Mantel$r2_Mantel <- round(df.Mantel$r2_Mantel, digit = 3)

  return(df.Mantel)
}

########## DONE! #










# /////////////////////////////////#
# phylogenetic_signal_pipeline ####
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# the functions below are used by the final function, phylogenetic_signal_pipeline, to calculate and plot the phylogenetic signal of a list of ps objects



# define a function that takes a nested list of ps objects ad returns an UMAP calculation for each list member
UMAP_on_ps_l_l <- function(ps_l_l) {
  # this function wil take a list of lists of phyloseq objects and return a list of lists of UMAP ordinations
  # ps_l_l = a list of lists of phyloseq objects. designed to take a phyloseq objectt were split by plant "Stress" condition


  # First transpose OTU table so that rows are samples and OTUs columns
  transposed_otu_table <- lapply(ps_l_l, function(x) {
    lapply(x, function(z) {

      # remove empty rows and columns
      output <- t(otu_table(z))
      output <- output[
        rowSums(output) > 0,
        colSums(output) > 0
      ]

      return(output)
    })
  })




  # Calculate dissimilarity indexes with bray curtis method
  bray_l <- lapply(transposed_otu_table, function(x) {
    lapply(x, function(z) {
      vegdist(z, method = "bray")
    })
  })

  # save it a a matrix
  bray_matrix_l <- lapply(bray_l, function(x) {
    lapply(x, function(z) {
      as.matrix(z)
    })
  })


  # Run UMAP on bray curtis dissimilarity matrix
  #PERO, MAKE THE NUMBER OF NEIGHTBORS = N-1 (N = NUMBER OF SAMPLES)
  set.seed(1234)
  UMAP_bray_l <- lapply(bray_matrix_l, function(x) {
    lapply(x, function(z) {
      umap(z, verbose = TRUE, n_epochs = 500)
    })
  })


  # output: the UMAP calculation result
  return(UMAP_bray_l)
}

# define function to extract umap coordinates and then merge it to the list of ps objects
add_UMAP_to_ps <- function(UMAP_on_ps_l_l_output, ps_l) {

  # this function wil take the output of the UMAP_on_ps_l_l function (a list of lists of UMAP ordinations) and add the coordinates to a ps object
  # UMAP_on_ps_l_l_output = the output of  UMAP_on_ps_l_l() function
  # ps_l = a single list a phyloseq objects. note that this is not a nested list: it's a lsit before it was split by "Stress"


  # Extract UMAP coordinates
  UMAP_coordinates_l_l <- lapply(UMAP_on_ps_l_l_output, function(x) {
    control_UMAP <- x$Control$layout
    MeJA_UMAP <- x$MeJA$layout
    SA_UMAP <- x$SA$layout

    # put control, meja and SA coordinates ona  single df to add to the ps object
    output <- as.data.frame(rbind(control_UMAP, MeJA_UMAP, SA_UMAP))

    colnames(output) <- c("UMAP_1", "UMAP_2")

    return(output)
  })

  # add umap coordiantes as metadata for the s object
  ps_with_umap_l <- mapply(function(x, y) {
    umap_coord_ps <- phyloseq::sample_data(y)

    x <- phyloseq::merge_phyloseq(x, umap_coord_ps)

    return(x)
  },
  x = ps_l,
  y = UMAP_coordinates_l_l
  )

  # output: ps_object with UMAP coordinates
  return(ps_with_umap_l)
}

#define function to calcualte PCoA and then add the PC1 and PC2 coordinates to the ps object
PCoA_calculate_extract_add_to_ps<-function(input_ps_l_l, original_ps_l){
  # input_ps_l_l = a nested list of phyloseq objects, that was split according "Stress
  # original_ps_l = a list of phyloseq objects, before they were split by "Stress"
  
  # calcualte PCoA
  set.seed(300)
  PCoA_stress_l_l <- lapply(input_ps_l_l, function (list_of_two)
    lapply(list_of_two, function (list_of_3)
      ordinate(list_of_3,
               method="PCoA",
               distance="bray",
               try=200,
               autotransform=TRUE)))
  
  
  
  #Now we extract the values from Axis1 and Axis2, returning a dataframe for each ps object
  PC1_PC2_df_l<-lapply(PCoA_stress_l_l, function(x){
    
    PC1_PC2_control <- as.data.frame(x$Control$vectors[,1:2])#ctrl
    PC1_PC2_MeJA <- as.data.frame(x$MeJA$vectors[,1:2]) #MeJA
    PC1_PC2_SA <- as.data.frame(x$SA$vectors[,1:2]) #SA
    
    PC1_PC2_df <- rbind(PC1_PC2_control, PC1_PC2_MeJA, PC1_PC2_SA)
    names(PC1_PC2_df)[1] <- "PCoA_1" #rename axis column
    names(PC1_PC2_df)[2] <- "PCoA_2" #rename axis column
    
    return(PC1_PC2_df)
    
  })
  
  
  
  
  # merge extracted coordinates with the rest of the phyloseq object
  ps_with_PCoA_coord<-mapply(function (x,y){
    
    PC1_PC2_ps <- phyloseq::sample_data(x)
    
    output<-phyloseq::merge_phyloseq(y, PC1_PC2_ps)
    
  },
  x = PC1_PC2_df_l,
  y = original_ps_l,
  SIMPLIFY = FALSE)
  
  
  
  return(ps_with_PCoA_coord)
  
  
  
  
  
}

# define function to calculate NMDS and then add the PC1 and PC2 coordinates to the ps object
NMDS_calculate_extract_add_to_ps<-function(input_ps_l_l, original_ps_l){
  
  # calculates NMDS on ensted list
  set.seed(300)
  NMDS_stress_BacFun_l <- lapply(input_ps_l_l, function (list_of_two)
    lapply(list_of_two, function (list_of_3)
      ordinate(list_of_3,
               method="NMDS",
               distance="bray",
               trymax=20,
               autotransform=TRUE)))
  
  
  
  
  
  #Now we extract the values from Axis1 and Axis2, returning a dataframe for each ps object
  PC1_PC2_df_l<-lapply(NMDS_stress_BacFun_l, function(x){
    
    PC1_PC2_control <- as.data.frame(x$Control$points)#ctrl
    PC1_PC2_MeJA <- as.data.frame(x$MeJA$points) #MeJA
    PC1_PC2_SA <- as.data.frame(x$SA$points) #SA
    
    PC1_PC2_df <- rbind(PC1_PC2_control, PC1_PC2_MeJA, PC1_PC2_SA)
    names(PC1_PC2_df)[1] <- "NMDS_1" #rename axis column
    names(PC1_PC2_df)[2] <- "NMDS_2" #rename axis column
    
    return(PC1_PC2_df)
    
  })
  
  
  
  
  # merge extracted coordinates with the rest of the phyloseq object
  ps_with_NMDS_coord<-mapply(function (x,y){
    
    PC1_PC2_ps <- phyloseq::sample_data(x)
    
    output<-phyloseq::merge_phyloseq(y, PC1_PC2_ps)
    
  },
  x = PC1_PC2_df_l,
  y = original_ps_l,
  SIMPLIFY = FALSE)
  
  
  
  return(ps_with_NMDS_coord)
  
}



# define function to calculate ASCA and then add the PC1 and PC2 coordinates to the ps object
ASCA_calculate_extract_add_to_ps<-function(input_ps_l_l, original_ps_l){
  
  # calculates ASCA on neseted list
  ASCA_stress_l_l <- lapply(input_ps_l_l, function (y)
    lapply(y, function (x){ 
      #put it into a format ASC will understand
      feature_table<-I(t(as.data.frame(x@otu_table)))
      asca_input<-as(sample_data(x), "data.frame")
      asca_input$ASVs<-feature_table
      
      # run asca
      mod<-asca(ASVs ~ Sp_abb_name + Block, data = asca_input) 
      
      return(mod)
      
    }
    )
  )
  
  
  
  
  #Now we extract the values from Axis1 and Axis2, returning a dataframe for each ps object
  PC1_PC2_df_l<-lapply(ASCA_stress_l_l, function(x){
    
    PC1_PC2_control <- as.data.frame(x$Control$projected$Sp_abb_name)[,1:2]#ctrl
    PC1_PC2_MeJA <- as.data.frame(x$MeJA$projected$Sp_abb_name)[,1:2] #MeJA
    PC1_PC2_SA <- as.data.frame(x$SA$projected$Sp_abb_name)[,1:2] #SA
    
    PC1_PC2_df <- rbind(PC1_PC2_control, PC1_PC2_MeJA, PC1_PC2_SA)
    names(PC1_PC2_df)[1] <- "ASCA_1" #rename axis column
    names(PC1_PC2_df)[2] <- "ASCA_2" #rename axis column
    
    return(PC1_PC2_df)
    
  })
  
  
  # print explained variance  
  explained_var<-lapply(ASCA_stress_l_l, function(x) lapply(x, function(y) y$explvar))
  explained_var_df_l<- lapply(explained_var, function(x)
    do.call("rbind", x))
  print("variance explained in ASCA")
  print(explained_var_df_l)
  
  
  
  # merge extracted coordinates with the rest of the phyloseq object
  ps_with_ASCA_coord<-mapply(function (x,y){
    
    PC1_PC2_ps <- phyloseq::sample_data(x)
    
    output<-phyloseq::merge_phyloseq(y, PC1_PC2_ps)
    
  },
  x = PC1_PC2_df_l,
  y = original_ps_l,
  SIMPLIFY = FALSE)
  
  
  
  return(ps_with_ASCA_coord)
  
}



# define a function to load and prepare the MY tree
load_prepare_MY_tree <- function(plot_Bra_tree = FALSE) {

  # 2.0b -load and prepare base substitution tree million year tree

  #Bra_tree <- read.nexus("./Plant_phylogeny/Results/final_tree.tre")
  #Bra_tree <- read.tree("./Plant_phylogeny/Supermatrix_Kasper_tree_renamed.tre") # new kaster tree
 # Bra_tree <- read.tree("./Plant_phylogeny/Supermatrix_Kasper_tree_renamed_rotated.tre") # new kaster tree
  Bra_tree <- read.tree("./17_Plant_phylogeny/Results/Supermatrix_Kasper_tree_renamed_rotated_rooted.tre") # new kaster tree
  
  # remove outgroup
  # Bra_tree <- drop.tip(Bra_tree, tip = "Aethionema_arabicum", trim.internal = TRUE, rooted = TRUE)

  # renaming tips
 # old_tips <- Bra_tree$tip.label
 # phylo_tips <- Bra_tree$tip.label
 # phylo_tips <- c(
#    "At", "Al", "Cas", "Tg", "Ec", "Mm", "Bv",
 #   "Lr", "Ds", "Bi", "Bn", "Sar", "Hi", "Br", "Bo_M", # attention on Bo_M or Bo.M;
#    "Dt", "Si", "It", "Co", "Lm", "Ia", "Es"
 # )

  # change for new names
 # Bra_tree$tip.label[match(old_tips, Bra_tree$tip.label)] <- phylo_tips

  # plot for cheking
  if (plot_Bra_tree == TRUE){
   plot(Bra_tree)
  }

  # output: MY brassicacee tree
  return(Bra_tree)
}

# define function to prepare a df used as input to ggplot
Abouheif_to_df <- function(Abouheif_and_moran_output, ps_to_df_output) {

  # Extract Obs and pvalues
  Abouheif_moran_Obs_Bac_l <- lapply_lapply(Abouheif_and_moran_output, function(x) x$obs)
  Abouheif_moran_pvalue_Bac_l <- lapply_lapply(Abouheif_and_moran_output, function(x) x$pvalue)

  # Convert the list to a named vector
  Abouheif_moran_obs_Bac_v <- unlist(Abouheif_moran_Obs_Bac_l, recursive = FALSE)
  Abouheif_moran_pvalue_Bac_v <- unlist(Abouheif_moran_pvalue_Bac_l, recursive = FALSE)

  # Stack the vector and create a data frame
  Abouheif_moran_obs_Bac_df <- as.data.frame(stack(Abouheif_moran_obs_Bac_v))
  Abouheif_moran_pvalue_Bac_df <- as.data.frame(stack(Abouheif_moran_pvalue_Bac_v))
  Abouheif_moran_Bac_df <- cbind(Abouheif_moran_obs_Bac_df, Abouheif_moran_pvalue_Bac_df) # merge by="ind" didn't work
  Abouheif_moran_Bac_df <- Abouheif_moran_Bac_df[-4] # remove repeated column

  # Rename the columns
  colnames(Abouheif_moran_Bac_df) <- c("value", "index", "pvalue")

  Abouheif_moran_Bac_df <- Abouheif_moran_Bac_df[, c("index", "value", "pvalue")] # re-order columns

  Abouheif_moran_Bac_df <- Abouheif_moran_Bac_df %>% separate(index, c("stress", "index")) # separate stress and test

  # Add column with traits
  trait <- colnames(ps_to_df_output$Control)[3:4] # only the last 2 column names
  trait <- rep(trait, length.out = dim(Abouheif_moran_Bac_df)[1])

  # add UMAP traits
  Abouheif_moran_Bac_df <- cbind(Abouheif_moran_Bac_df, trait)

  # Remove stress trait rows
  Abouheif_moran_Bac_df <- subset(Abouheif_moran_Bac_df, trait != "Stress")

  # Keep only traits mean (remove SD)
  Abouheif_moran_Bac_df <- Abouheif_moran_Bac_df[!grepl("_sd", Abouheif_moran_Bac_df$trait), ]

  # Modify index names
  Abouheif_moran_Bac_df <- Abouheif_moran_Bac_df %>% mutate(index = str_replace(index, "moran", "Moran's I"))
  Abouheif_moran_Bac_df <- Abouheif_moran_Bac_df %>% mutate(index = str_replace(index, "abouheif", "Abouheif's Cmean"))

  return(Abouheif_moran_Bac_df)
}

# function to preapre a df for the plot_traits_on_tree() function
prepare_to_table.phylo4d <- function(ps_to_df_output, p_adjusted_df, significant_signal_only) {


  # this functions will wrangle the data from the output of ps_to_df to an inpute for table.phylo4d, filtering to the traits that were significant on p-adjusted abuheif test or not

  # ps_to_df_output = output of the ps_to_df() function
  # p_adjusted_df = df with adjusted p values of the phylogenetic signal test
  # significant_signal_only = Logical (TRUE of FALSE) indicating if you want to keep only the trats that were significant in the phylogenetic test


  # significant groups and stresses:
  sig_groups_abouf <- p_adjusted_df[p_adjusted_df$pvalue < 0.05, ][, c(1, 2, 6)]

  # significant traits
  sig_traits_abouhf <- paste(sig_groups_abouf$trait, sig_groups_abouf$stress, sep = "_")





  # if you asked to keep only the traits that were significant on the abuheif test, return only the abuheif test
  if (significant_signal_only == TRUE) {
    # a much smaller table; only includes taxons that were significant
    ps_to_df_output_subset <- ps_to_df_output[names(ps_to_df_output) %in% sig_groups_abouf[, 1]]
  } else {
    # the full dataset
    ps_to_df_output_subset <- ps_to_df_output
  }






  output1 <- lapply(ps_to_df_output_subset, function(z) {
    # this function will change some column names to reflext associates stresses
    # it takes the output of ps_to_df and then...
    # removes "mean"from trait names, adding the "stress"instead
    lapply(z, function(y) {
      new_col_names <- paste(colnames(y),
        y[[1, 2]], # the stress in the df list)
        sep = "_"
      )

      new_col_names2 <- gsub(pattern = "_mean", replacement = "_temp", x = new_col_names)

      colnames(y)[3:4] <- new_col_names2[3:4]

      return(y)
    })
  })








  # unlist, remove uncessary metadata columns that were replicated
  output2 <- lapply(output1, function(x) do.call(cbind, x)[, c(1, 3:4, 7:8, 11:12)])


  # adjust variable names to match names of lists
  output3 <- mapply(function(z, y) {
    colnames(z) <- as.factor(gsub(
      pattern = "temp",
      replacement = y,
      x = colnames(z)
    ))

    colnames(z)[1] <- "Sp_abb_name"
    return(z)
  },
  z = output2,
  y = names(output2),
  SIMPLIFY = FALSE
  )


  # join all list of dfs into a single df
  output4 <- plyr::join_all(output3, by = "Sp_abb_name")

  # only keep the column names (trats) that had a significant abouheif results
  output5 <- dplyr::select(output4, Sp_abb_name, sig_traits_abouhf)

  # if you asked to keep only the traits that were significant on the abuheif test, return only the abuheif test
  if (significant_signal_only == TRUE) {
    output6 <- output5
  } else {
    output6 <- output4
  }


  return(output6)
}

# function to plot trait next to phylogeny
plot_traits_on_tree <- function(df_to_table.phylo4d_output, phylogenetic_tree) {
  # This funtion will take a df in the format of the ps_to_df() fucntion and match such emtadata agasint the phylogenetic tree
  # it has two arguments: the output of ps_to_df() and a phylogenetic tree like Bra_tree
  # the output are the metrics for the phylogenetic signal (Abouheif and abouheif )

  # prepare data
  test <- as.data.frame(df_to_table.phylo4d_output, row.names = NULL) # change to df, otherwise weird errors
  row.names(test) <- test[, 1] # add rownames
  test <- test[, -1] # remove 'Sp' column

  # make adephylo object
  phylotraits <- phylo4d(phylogenetic_tree, test) # from phylobase package

  table.phylo4d(phylotraits)

  output <- recordPlot()

  return(output)
}


# runs a phylogenetic principal component on the traits
ppca_traits_on_tree <- function(df_to_table.phylo4d_output, phylogenetic_tree) {
  # This funtion will take a df in the format of the ps_to_df() fucntion and match such emtadata agasint the phylogenetic tree
  # it has two arguments: the output of ps_to_df() and a phylogenetic tree like Bra_tree
  # the output are the metrics for the phylogenetic signal (Abouheif and abouheif )
  
  # prepare data
  test <- as.data.frame(df_to_table.phylo4d_output, row.names = NULL) # change to df, otherwise weird errors
  row.names(test) <- test[, 1] # add rownames
  test <- test[, -1] # remove 'Sp' column
  
  # make adephylo object
  phylotraits <- phylo4d(phylogenetic_tree, test) # from phylobase package
  
  #save "phylotraits" in the global enviroment so we can make some plots
  assign("phylotraits", phylotraits, envir=globalenv())
  
  #calcualte ppca
  ppca_flattened <- ppca(phylotraits, scannf=FALSE, nfposi = 2, nfnega =0, method="oriAbouheif")
  

  #save pltos and ouput  
  plot(ppca_flattened, useLag=TRUE)
  output1 <- recordPlot()
  
  print(ppca_flattened)
  summary(ppca_flattened)
  
  output2<-ppca_flattened
  output3<-summary(ppca_flattened)
  
  #define output
  output <-list("ppca_plot" = output1,
               "ppca_obj" = output2,
               "ppca_summary" = output3)
  
  # remove phylotraits from global eviroment to aovid confusion
  rm(phylotraits, envir=globalenv())
  
  return(output)
}


#define a function that will take a single df from random_ASV_phyloSig2 and add a new column that separates the higer PC from the lower PC
tag_higer_lower_PC<-function(random_ASV_phyloSig2_output_df){
  # random_ASV_phyloSig2_output_df = one single df from random_ASV_phyloSig2()
  # returns the same df, with a new column that discriminates PC1 and PC2 as higehr or lower in phylogenetic signal for that ordination
  
  #extract vectors
  pc1<-dplyr::filter(random_ASV_phyloSig2_output_df, grepl("_1_", trait))$value
  pc2<-dplyr::filter(random_ASV_phyloSig2_output_df, grepl("_2_", trait))$value
  
  
  #define true/false
  comp1<-pc1>=pc2
  comp2<-pc1<pc2
  PC_comparison<-c(rbind(comp1,comp2))
  
  
  #create new column
  random_ASV_phyloSig2_output_df_mutated<-mutate(.data = random_ASV_phyloSig2_output_df, 
                                                 higher_lower_PC = ifelse(test =  PC_comparison == TRUE, yes = "Higher_PC", no = "Lower_PC" ))
  
  return(random_ASV_phyloSig2_output_df_mutated)
  
}



# define function to run whole pipelie, consolidating several custom functions
phylogenetic_signal_pipeline <- function(ps_l, 
                                         phylogenetic_test = "Abouheif's Cmean", 
                                         ordination_method = "PCoA", 
                                         significant_signal_only = FALSE,
                                         return_plot = TRUE,
                                         return_lambda = FALSE,
                                         higher_PC_only = FALSE,
                                         show_abouheif_histogram = FALSE,
                                         ...) {
  
  # These are the arguments of the phylogenetic_signal_pipeline function:
    # ps_l = a  list of phyloseq objects. list elements must have a name, which will be pasted after the UMAP dimensions on the plot
    # phylogenetic_test = a quoted charatehr variable, allowed values are "Moran's I" or "Abouheif's Cmean"
    # significant_signal_only = logical, should you keep in the plots only the traits that had a significant test value? default is FALSE
    # ordination_method = srng, which ordination method to calculat the centroid values from. allowed values are "PCoA", "NMDS", "UMAP" and "ASCA". default is "PCoA"
    # return_plot = should the function return a plot with hylosig values values? default is TRUE 
    # higher_PC_only = logical, should the plot retun only the higest of the two PCs?  default is FALSE. 
    # show_abouheif_histogram = logical, should you show a histogram of the abouheif metric against the random pairs? default is FALSE
  
  # This function returns:
    # plot_abouheif = a barplot of abouheif value under each stress and each phyloseq object
    # trait_on_tree =  a trait_on_tree plot, showing the ordinaton coordinates as black/white circles in the phylogenetic tree
    # ppca_on_tree = a ppca of all stress conditions and phyloseq objects
    # df_phylosig = a df with phylogenetic signal values for each ordination axis and phyloseq object
    # ordination_coordinates_l = a phyloseq object with the input data, also including the coordinates of each sample in the ordination
  
  
  # load brssica MY tree
  Bra_tree <- load_prepare_MY_tree()

  # turn a ps_l into a nested list according variable "Strees"
  ps_l_l <- lapply(ps_l, function(x) {
    phyloseq_sep_variable(x, variable = "Stress")
  })

  

  # running ordination_method = "UMAP"
  if (ordination_method == "UMAP") {
    
    # calculate umap
    Umap_output_l_l <- UMAP_on_ps_l_l(ps_l_l)
    
    # add umap coordinates to ps_l
    ps_with_coordinates <- add_UMAP_to_ps(UMAP_on_ps_l_l_output = Umap_output_l_l, ps_l = ps_l)
    
    # define metadata columns to keep
    metadata_columns<-c("Sp_abb_name",
                        "Stress",
                        "UMAP_1",
                        "UMAP_2")
  }
  
  
  
  
  # running ordination_method = "PCoA"
  if (ordination_method == "PCoA") {

    #calcualte PCoA and put PCoA_1 and PCoA_2 coordinates on ps object
    ps_with_coordinates<-PCoA_calculate_extract_add_to_ps(input_ps_l_l = ps_l_l, original_ps_l = ps_l)
    
    # define metadata columns to keep
    metadata_columns<-c("Sp_abb_name",
                        "Stress",
                        "PCoA_1",
                        "PCoA_2")
  
  }
  
  
  
  
  
  # running ordination_method = "NMDS"
  if (ordination_method == "NMDS") {
    
    #calcualte PCoA and put PCoA_1 and PCoA_2 coordinates on ps object
    ps_with_coordinates<-NMDS_calculate_extract_add_to_ps(input_ps_l_l = ps_l_l, original_ps_l = ps_l)
    
    # define metadata columns to keep
    metadata_columns<-c("Sp_abb_name",
                        "Stress",
                        "NMDS_1",
                        "NMDS_2")
    
  }
  
  
  
  
  # running ordination_method = "ASCA"
  if (ordination_method == "ASCA") {
    
    #calcualte PCoA and put PCoA_1 and PCoA_2 coordinates on ps object
    ps_with_coordinates<-ASCA_calculate_extract_add_to_ps(input_ps_l_l = ps_l_l, original_ps_l = ps_l)
    
    # define metadata columns to keep
    metadata_columns<-c("Sp_abb_name",
                        "Stress",
                        "ASCA_1",
                        "ASCA_2")
    
  }
  
  
  # turn ps_obj into a df for phylogenetic test
  df_plant_l <- lapply(ps_with_coordinates, function(x) {
    ps_to_df(x, colstokeep = metadata_columns)
  })

  
  
  
  
  # run the Abouheif_and_moran tests of metadata on plat phylogenies
  Abouheif_and_moran <- lapply(df_plant_l, function(y) {
    lapply(y, function(x) {
      plantdf_to_Abouheif_and_moran(
        ps_to_df_output = x,
        phylogenetic_tree = Bra_tree,
        show_abouheif_histogram = show_abouheif_histogram
      )
    })
  })

  
  # execute Abouheif_to_df fucntion over lsit
  df_to_ggplot_l <- mapply(function(y, x) {
    Abouheif_to_df(
      Abouheif_and_moran_output = x,
      ps_to_df_output = y 
    )
  },
  x = Abouheif_and_moran,
  y = df_plant_l,
  SIMPLIFY = FALSE
  )







  # adjust variable names to match names of lists
  df_to_ggplot_l <- mapply(function(z, y) {
    z$trait <- as.factor(gsub(
      pattern = "mean",
      replacement = y,
      x = z$trait
    ))
    return(z)
  },
  z = df_to_ggplot_l,
  y = names(df_to_ggplot_l),
  SIMPLIFY = FALSE
  )


  # put all list slices into a single df
  df_to_ggplot <- bind_rows(df_to_ggplot_l, .id = "column_label")

  # remove moran's tests
  df_to_ggplot <- filter(df_to_ggplot, index == phylogenetic_test)

  # add a new column to the df, differentiating the higher PC from the lwoer PC
  df_to_ggplot<-tag_higer_lower_PC(df_to_ggplot)
  
  
  
# run the lambda and K tests of metadata on plat phylogenies, PC_1
  lambda_and_k_dim1 <- lapply(df_plant_l, function(y) {
    lapply(y, function(x) {
      phylogenetic_lambda_K(
        ps_to_df_output = x,
        selected_trait = colnames(x)[3],
        phylogenetic_tree = Bra_tree
      )
    })
  })
  
  
  
# run the lambda and K tests of metadata on plat phylogenies, PC_2  
  lambda_and_k_dim2 <- lapply(df_plant_l, function(y) {
    lapply(y, function(x) {
      phylogenetic_lambda_K(
        ps_to_df_output = x,
        selected_trait = colnames(x)[4],
        phylogenetic_tree = Bra_tree
      )
    })
  })
  
  
  
  #extract lambda values into a df_l
  lambda_df_l<-mapply(function(x,y){
    
    
    # get results into a small df
    lambda_dim1<-lambda_k_results(input_l = x)
    lambda_dim2<-lambda_k_results(input_l = y)
    
    # label dimensions
    lambda_dim1$trait<-"Dimension_1_"
    lambda_dim2$trait<-"Dimension_2_"
    
    # intercaalte dimension 1 and 2, like on abuheif tests
    # note: here we discard Blomberg's K
    lambk_df<-rbind(lambda_dim1[1,], # could not automate in time
                    lambda_dim2[1,],
                    lambda_dim1[2,],
                    lambda_dim2[2,],
                    lambda_dim1[3,],
                    lambda_dim2[3,])
    
    
    return(lambk_df) 
  },
  x = lambda_and_k_dim1,
  y = lambda_and_k_dim2,
  SIMPLIFY = FALSE)
  
  #consolidate df_l into a single df
  lambda_df<-do.call("rbind", lambda_df_l)
  
  # tag PCs 1 and 2 as higher or lower
  lambda_df<-tag_higer_lower_PC(lambda_df)
  
  # if you want the labda test values, overwrite the abouheif.moran's index
  if (return_lambda == TRUE){
    
    index<-"Pagel's lambda"

    df_to_ggplot$index<-lambda_df$index
    df_to_ggplot$value<-lambda_df$value
    df_to_ggplot$pvalue<-lambda_df$pvalue
    df_to_ggplot$higher_lower_PC<-lambda_df$higher_lower_PC
  } 
  
  
  
  if (higher_PC_only == TRUE){
    df_to_ggplot<- filter(df_to_ggplot, higher_lower_PC == "Higher_PC")
    
    
  } 
  
  
  
  
  # adjsut p values with fdr
  df_to_ggplot$pvalue_adjust <- p.adjust(p = df_to_ggplot$pvalue, method = "fdr")
  
  
  
  
  # filter the DF if higher_PC is TRUE
  if (higher_PC_only == TRUE){
    df_to_ggplot<- filter(df_to_ggplot, higher_lower_PC == "Higher_PC")
    
    
  } 


  # Set axis looks
  theme_set(theme_bw())
  axis_looks <- theme(axis.text.x = element_text(
    colour = "black", size = 10,
    face = "bold"
  )) +
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold")) +
    theme(axis.title = element_text(size = 15, face = "bold"))


  Color <- brewer.pal(4, "BrBG")
  

  
  # if plots are requested, execute the ploting functions, other wide jsut return df
  if (return_plot == TRUE){

  # generate plot
  p.summary <- ggplot(
    data = df_to_ggplot,
    mapping = aes(x = trait, y = value, fill = index)
  ) +
    geom_bar(
      stat = "identity",
      position = "dodge"
    ) +
    scale_fill_manual(values = c("#80CDC1", "#018571")) +
    geom_hline(aes(yintercept = 0), linetype = "solid") +
    coord_flip() +
    theme(axis.text.y = element_text(size = 6)) +
    scale_x_discrete(limits = rev) +
    facet_wrap(~stress) +
    geom_text(aes(label = ifelse(test = pvalue < 0.05, yes = pvalue, no = ""))
    ) +
    geom_text(aes(label = ifelse(test = pvalue_adjust < 0.05, yes = "*", no = "")),
              hjust = 0.5, size = 11, inherit.aes = TRUE) +
    axis_looks
  # print plot
  print(p.summary)

  # prepare df for table.phylo4d
  phylo4d_df <- prepare_to_table.phylo4d(
    ps_to_df_output = df_plant_l,
    p_adjusted_df = df_to_ggplot,
    significant_signal_only = significant_signal_only
  )

  #check df
  print(phylo4d_df)

  if (length(ps_l) <30){
  
  # plot table.phylo4d
  trait_on_tree <- plot_traits_on_tree(
    df_to_table.phylo4d_output = phylo4d_df,
    phylogenetic_tree = Bra_tree
  )
  
  } else {
    trait_on_tree<-"warning: more than 30 ps objects cannot be labeled in a trait tree with this function"  
  }
  
  # separate ordinations of control meja ans SA for separate PPCAs in Control MeJA and SA
  ppca_input_l<-list(
  "ppca_input_Control" = dplyr::select(phylo4d_df, Sp_abb_name, ends_with("Control")),
  "ppca_input_MeJA" = dplyr::select(phylo4d_df, Sp_abb_name, ends_with("MeJA")),
  "ppca_input_SA" = dplyr::select(phylo4d_df, Sp_abb_name, ends_with("SA")))
  
  # run ppca
#  ppca_on_tree <- ppca_traits_on_tree(
#    df_to_table.phylo4d_output = phylo4d_df,
#    phylogenetic_tree = Bra_tree
#  )
  
  # run ppca_l
  ppca_on_tree <- lapply(ppca_input_l, function (x)
    ppca_traits_on_tree(
    df_to_table.phylo4d_output = x,
    phylogenetic_tree = Bra_tree
  ))

  # define output
  output <- list(
    "plot_abouheif" = p.summary,
    "trait_on_tree" = trait_on_tree,
    "ppca_on_tree" = ppca_on_tree,
    "df_phylosig" = df_to_ggplot,
    "ordination_coordinates_l" = ps_with_coordinates
  )
  

} else {
  output <-  list( "df_phylosig" = df_to_ggplot,
                   "ordination_coordinates_l" = ps_with_coordinates)
}

  return(output)
}





# /////////////////////////////////#
# Permanova and PCoA-based filtering ####
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#

# permanova and beta dispersion tests for  bray_curts ~ Sp_full_name + Stress + Sp_full_name*Stress + Block
permanova_per_taxonomic_group<-function(ps_obj_l){
  # ps_obj_l = a list of phyloseq objects, split by taxonomic group
  # output: A list of permanova results ( ~ Sp_full_name + Stress + Sp_full_name*Stress + Block) and...
  #         A list of beta dispersion results (Species, BLock, Stress)
  # runs quickly on the HPC with  parallel::mclapply (many cores, low memory)
  
  
  
  ps_obj_l<- lapply(ps_obj_l, function(x){
    
    x@sam_data$Sp_Lineage_Walden_extended<-
      if_else(condition = x@sam_data$Sp_abb_name %in% c("Ia", "Lm", "Co", "Bi"), 
              true ="lineage_II_extended",
              false =  x@sam_data$Sp_Lineage_Walden)
    
    return(x)
    
  })
  
  
  #calculate permanova at individual stress conditions
  set.seed(303848)
  permanova_l <- lapply(ps_obj_l, function (x)
    adonis2(formula = phyloseq::distance(t(otu_table(x)), method="bray") 
            ~ Sp_full_name + Stress + Sp_full_name*Stress + Block,
            data = as(sample_data(x),"data.frame"), # changing with as.data.frame is insufficient
            permutations = 999, by = "terms"))# reduced to facilidate computation
  
  
  beta_disper_l <- lapply(ps_obj_l, function (x){
    
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
  })
  
  
  
  output<-list("permanova_l" = permanova_l,
               "beta_disper_l" = beta_disper_l)
  
  return(output)
  
}

# this function will extract r and R values
permanova_output_to_plot_input2<-function(permanova_output_l){
  #permanova_output_l = the list of permanovas ("$permanova_l") from the function permanova_per_taxonomic_group()
  # the output is a df with p and R2 values for block, species, sptress and intereaction
  
  
  
  p_values <-purrr::map(permanova_output_l, 5)
  R2_values <- purrr::map(permanova_output_l, 3)
  list_lenght<-length(p_values)
  
  #make vectors with data of interest
  vector_of_p_stress <- as.vector(unlist(purrr::map(p_values,1)))
  vector_of_p_block <- as.vector(unlist(purrr::map(p_values,2)))
  
  vector_of_R2_Sp <- as.vector(unlist(purrr::map(R2_values,1)))
  vector_of_R2_Stress <- as.vector(unlist(purrr::map(R2_values,2)))
  vector_of_R2_Block <- as.vector(unlist(purrr::map(R2_values,3)))
  vector_of_R2_Sp_x_Stress <- as.vector(unlist(purrr::map(R2_values,4)))
  vector_of_R2_residual <- as.vector(unlist(purrr::map(R2_values,5)))
  
  
  vector_of_p_Sp <- as.vector(unlist(purrr::map(p_values,1)))
  vector_of_p_Stress <- as.vector(unlist(purrr::map(p_values,2)))
  vector_of_p_Block <- as.vector(unlist(purrr::map(p_values,3)))
  vector_of_p_Sp_x_Stress <- as.vector(unlist(purrr::map(p_values,4)))
  vector_of_p_residual <- as.vector(unlist(purrr::map(p_values,5)))
  
  
  # make a vector with all p values
  vector_p <- c(vector_of_p_Sp,
                vector_of_p_Stress,
                vector_of_p_Block,
                vector_of_p_Sp_x_Stress)
  #  vector_of_p_residual)
  
  
  vector_R2 <- c(vector_of_R2_Sp,
                 vector_of_R2_Stress,
                 vector_of_R2_Block,
                 vector_of_R2_Sp_x_Stress)
  #  vector_of_R2_residual)
  
  #defines factors
  factor_text<-c(rep("Sp_full_name", list_lenght),
                 rep("Stress", list_lenght),
                 rep("Block", list_lenght),
                 rep("Sp_x_Stress", list_lenght))
  #   rep("residual", list_lenght))
  
  
  species_name<-as.factor(names(c(R2_values)))
  
  
  #join vectors in a dataframe
  
  p_df <-  data.frame(vector_p, # p values for the stress and then for the block
                      species_name, # names of species in list
                      as.factor(factor_text)) # Stress and block as factors
  
  R2_df <-  data.frame(vector_R2, # p values for the stress and then for the block
                       species_name, # names of species in list
                       as.factor(factor_text)) # Stress and block as factors
  
  output_df <- merge(p_df, R2_df)
  
  #change anmes of columns in dataframe
  names(output_df)<-c("taxonomic_group","factor", "pvalue", "R2")
  
  #output
  return(output_df)
  
}

#this function will extract the varaince each PCoA explained, and also plot 500 PCoAs
explained_PCoA_Axis_per_taxonomic_group<-function(ps_obj_l){
  # ps_obj_l = a list of phyloseq objects, split by taxonomic group
  # the outut is: a list of PCoA ordiantions (one per phyloseq object) and...
  #               A list of PCoA plots (one per phyloseq object) and...
  #               A list of explained variances (one per phyloseq object)
  
  # add extended lineage as metadata
  ps_obj_l<- lapply(ps_obj_l, function(x){
    
    x@sam_data$Sp_Lineage_Walden_extended<-
      if_else(condition = x@sam_data$Sp_abb_name %in% c("Ia", "Lm", "Co", "Bi"), 
              true ="lineage_II_extended",
              false =  x@sam_data$Sp_Lineage_Walden)
    
    return(x)
    
  })
  
  
  
  
  PCoA_l<- lapply(ps_obj_l, function(x)
    ordinate(physeq = x,
             method="PCoA",
             distance="bray",
             autotransform=TRUE))
  
  
  ploted_ordination<- mapply(function (x,y)
    plot_ordination(physeq = x,
                    ordination = y, 
                    color = "Sp_Lineage_Walden_extended",
                    shape = "Stress")+
      scale_shape_manual(values = c(19,23,3)),
    x = ps_obj_l,
    y = PCoA_l,
    SIMPLIFY = FALSE)
  
  explained_variance<-lapply(ploted_ordination, function(z)
    c(z$labels$x, z$labels$y))
  
  output<-list("PCoA_l" = PCoA_l,
               "ploted_ordination" = ploted_ordination,
               "explained_variance" = explained_variance)
  
  return(output)
  
  
}

#this functions puts the % explaned variance of a PCoA intoa df
PCoA_plot_l_to_pct_variance_df<-function(explained_PCoA_output_l){
  # explained_PCoA_output_l = the output of explained_PCoA_Axis_per_taxonomic_group
  # the otput is the % varaince explained of a PCoA, in a df
  
  axis_pct_tax_group<-as.data.frame(names(explained_PCoA_output_l$explained_variance))
  names(axis_pct_tax_group)<-"taxonomic_group"
  
  axis_pct_tax_group$PCoA_1<-lapply(explained_PCoA_output_l$explained_variance, function(x) x[1])%>%unlist
  axis_pct_tax_group$PCoA_1<-gsub(x = axis_pct_tax_group$PCoA_1, pattern = "Axis.1   \\[", replacement = "")
  axis_pct_tax_group$PCoA_1<-as.numeric(gsub(x = axis_pct_tax_group$PCoA_1, pattern = "\\%]", replacement = ""))
  
  axis_pct_tax_group$PCoA_2<-lapply(explained_PCoA_output_l$explained_variance, function(x) x[2])%>%unlist
  axis_pct_tax_group$PCoA_2<-gsub(x = axis_pct_tax_group$PCoA_2, pattern = "Axis.2   \\[", replacement = "")
  axis_pct_tax_group$PCoA_2<-as.numeric(gsub(x = axis_pct_tax_group$PCoA_2, pattern = "\\%]", replacement = ""))
  
  return(axis_pct_tax_group)
  
}


# scree plot of PCoA variance, with broke stick (ev = PCoA eigenvalue), sourced in online forum
evplot = function(ev) {  
  # Broken stick model (MacArthur 1957)  
  n = length(ev)  
  bsm = data.frame(j=seq(1:n), p=0)  
  bsm$p[1] = 1/n  
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))  
  bsm$p = 100*bsm$p/n  
  # Plot eigenvalues and % of variation for each axis  
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))  
  barplot(ev, main="Eigenvalues", col="bisque", las=2)  
  abline(h=mean(ev), col="red")  
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")  
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE,   
          main="% variation", col=c("bisque",2), las=2)  
  legend("topright", c("% eigenvalue", "Broken stick model"),   
         pch=15, col=c("bisque",2), bty="n")  
  par(op)  
} 