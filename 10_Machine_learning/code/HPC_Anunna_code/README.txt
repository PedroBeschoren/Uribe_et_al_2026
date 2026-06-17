README - hpc_aNUNNA_CODE: Running Machine leanining R scripts on the HPC

to run these scripts on the HPC, they must all be on the Rpoject folder
run the .sb files with sbatch. the sb script will call rscript in the HPC after loading the R module, which will let you run the code insite the .R files

boruta_8.R:script to define ASVs (bacterial and fungal) that separate control samples from MeJA samples, and control samples from SA samples. it is run species by species in forests wtih 400k trees
CV_boruta400k_replicated25_MeJA_33sp.R:script that take the important RF features that separate control and MeJA to calculate precision, accuracy, etc. this is performed 25 times to account for random sampling when splitting training and testing models. This script also train the model based on boruta-selected features
CV_boruta400k_replicated25_SA_33sp.R:script that take the important RF features that separate control and SA to calculate precision, accuracy, etc. this is performed 25 times to account for random sampling when splitting training and testing models. This script also train the model based on boruta-selected features

boruta_8.sb:  runs boruta_8.R in the HPC Anunna at WUR (after cloning github repository, installing packages)
CV_boruta400k_replicated25_MeJA_33sp.sb: runs CV_boruta400k_replicated25_MeJA_33sp.R in the HPC Anunna at WUR (after cloning github repository, installing packages)
CV_boruta400k_replicated25_SA_33sp.sb: runs CV_boruta400k_replicated25_SA_33sp.R in the HPC Anunna at WUR (after cloning github repository, installing packages)