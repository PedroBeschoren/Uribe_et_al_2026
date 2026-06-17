# This script lists all libraries that sould be loaded, and then loads them all 

# this will install any packages present in the in the lockfile (renv.lock) but missing in your library
# source("./renv/restore_with_renv.R")


# trying to add this at the beginning of the script so Renv gets activated and appears in 'Packages' 
renv::load()

# lists all base libraries we need. li
packages <- c("devtools",# needed to install some packages
              "BiocManager", # needed to install some packages
              "remotes",# needed to install some packages
              "ggplot2",
              "dplyr",
              "tibble", 
              "rbiom", #needed to load large 16S biom files
              "stringr",# to wrangle string vectors
              "Boruta", # for random forest feature selection
              "mlbench", # for random forest
              "caret", # for random forest
              "randomForest", # for random forest
              "tidyr",
              "vegan", # for several essential statistical tests
              "forcats",
              "ggrepel", # to avoid legends overlapping in your plot
              "ggpubr",
              "igraph", # calculates entowrk metrics and manipulates netowrk objects
              "metagMisc", #  lets you create lists of split phyloseq objects
              "pheatmap", # heatmaps for deseq2
              "agricolae",# includes some anova post-hoc options
              "minpack.lm", #lets you do som HSD tests, output is a nive table
              "Hmisc", # for neutral models
              "spaa", # need to install Ecoutils
              "stats4",# for neutral models
              "car", #for levene test
              "metacoder", # plots heat trees
              "purrr", # has map() to select table elements
              "viridis", # prety colors
              "phyloseq", # essential to produce phyloseq objects with OTU, metadata, and taxa info in one single place
              "microbiome", # has convinient data wrangling functions
              "decontam", # to make use of no-template blank DNA extractions
              "metagenomeSeq", # to normalize library sizes without rarefying them
              "DESeq2", # used for differential abundance analysis
              "WGCNA", # needed for eigen_correlation (), allowing you to correlate metadata to network modules
              "pairwiseAdonis",
              "tyRa",  #need for neutral models
              "metagMisc",  #n lets you split a phyloseq object in a list
              "EcolUtils",   #pairwise adonis fucntion
              "MicEco",  #for venn diagrams on phylosseq  objects
              "SpiecEasi") # builds the sparse networks


# set the path of the libraries for the cluster
.libPaths("./renv/library/R-4.1/x86_64-w64-mingw32")
#.libPaths("C:/RProjects/Family_experiment/renv/library/R-4.1/x86_64-w64-mingw32") use this for MeJa pilot


## Now load these libraries 
lapply(
  packages,
  FUN = function(x) {
    require(x, character.only = TRUE)
  }
)

# remove package list to avoid visual pollution of your work environment
rm(packages)

#load some relevant custom functionseveyrone uses, like lapply_lapply
source("./renv/functions/lapply_lapply.R")

#let's defined the colour we will be using for our plots

colour_v <- c("#F0E442", "#88CCEE", "#E69F00", "#009E73", "#999999") #Bac #Fun #MeJA #SA #Control
