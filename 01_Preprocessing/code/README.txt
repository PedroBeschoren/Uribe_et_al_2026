script Files in this folder:

1_dada2fix_loading_decontmination_normalization.Rmd: key script to load raw data, process it, and generate filtered phyloseq objects
functions_loading_and_decontamination.R: custom function used to load, filter and decontaminate the data (some fucntions are not used)
host_DNA_contamination.R: quick checks on plastids, mithochondrial, and plant/protozoan ITS contaminations. does not have to be run.
merging_metadata_files.R: generates the mapping files for 16S and ITS. does not have to be run unless updating the mapping file
DNA_buffers_and_yield.R: calculations on the amount of buffer necessary for DNA extraction, according root weight. does not ahve to be run