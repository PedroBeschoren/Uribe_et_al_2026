#!/bin/bash -login
 
#SBATCH --time=00:01:00				  	### limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1					      ### number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=16			  ### number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=32G					        ### memory required per node - amount of memory (in bytes)
#SBATCH --job-name="tree"			### you can give your job a name for easier identification (same as -J)
#SBATCH --output="tree_%j.output"


# activates qiime2 environment
conda activate qiime2-2022.8

#From Qiime2(https://docs.qiime2.org/2022.8/tutorials/phylogeny/): 
#This pipeline will start by creating a sequence alignment using MAFFT,
#after which any alignment columns that are phylogenetically uninformative or ambiguously aligned will be removed (masked).
#The resulting masked alignment will be used to infer a phylogenetic tree and then subsequently rooted at its midpoint.
#Output files from each step of the pipeline will be saved.
#This includes both the unmasked and masked MAFFT alignment from q2-alignment methods, and both the rooted and unrooted phylogenies from q2-phylogeny methods.


#qiime phylogeny align-to-tree-mafft-fasttree \
#  --i-sequences rep-seqs.qza \
#  --output-dir mafft-fasttree-output

#qiime tools export \
#  /mafft-fasttree-output/unrooted-tree.qza \
#  --output-dir exported-tree

# delete this later:
# imports representative sequences in FASTA format
#qiime tools import \
#  --input-path /lustre/nobackup/INDIVIDUAL/arago004/tree_Family_filteredBac_08-11/refseq_filt_Bac.fna \
#  --output-path /lustre/nobackup/INDIVIDUAL/arago004/tree_Family_filteredBac_08-11/refseq_filt_Bac.qza \
#  --type 'FeatureData[Sequence]'

# pipeline to create phylogenetic tree
#qiime phylogeny align-to-tree-mafft-fasttree \
#  --i-sequences /lustre/nobackup/INDIVIDUAL/arago004/tree_Family_filteredBac_08-11/refseq_filt_Bac.qza \
#  --o-alignment /lustre/nobackup/INDIVIDUAL/arago004/tree_Family_filteredBac_08-11/refseq_filt_Bac_aligned.qza \
#  --o-masked-alignment /lustre/nobackup/INDIVIDUAL/arago004/tree_Family_filteredBac_08-11/refseq_filt_Bac_aligned_masked.qza \
#  --o-tree /lustre/nobackup/INDIVIDUAL/arago004/tree_Family_filteredBac_08-11/refseq_filt_Bac_tree_unrooted.qza \
#  --o-rooted-tree /lustre/nobackup/INDIVIDUAL/arago004/tree_Family_filteredBac_08-11/refseq_filt_Bac_tree_rooted.qza \
#  --verbose
  
# export tree from .qza to .nwk format 
qiime tools export \
  --input-path /lustre/nobackup/INDIVIDUAL/arago004/tree_Family_filteredBac_08-11/refseq_filt_Bac_tree_rooted.qza \
  --output-dir /lustre/nobackup/INDIVIDUAL/arago004/tree_Family_filteredBac_08-11/refseq_filt_Bac_tree_rooted 
  
conda deactivate


