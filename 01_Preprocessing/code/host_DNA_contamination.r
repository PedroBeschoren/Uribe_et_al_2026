#quantify host DNA contamination

# loads the bulk of the libraries necessary for the family experiment project
source("./src/Amplicon_sequencing/Functions/0_load_libraries.r")
source("./src/Amplicon_sequencing/Functions/loading_and_decontamination_custom_functions.R")



######## ----- 16S data ---------- ########

#load 16s data from biome file and makke a phyloseq object (requires ~ 8gig ram free)

physeq_16S<-load_16S_biom("./data/Amplicon_sequencing/16S_phyloseq_input_FeatureTable_metadata_taxonomy.biom")



# removes the horrible ASV naming from qiime2 into ASV numbers
physeq_all<-backup_and_rename(physeq_all) 






### a quick check to detect the presence of  o__Chloroplast or f__Mitochondria

# set some vectors with samples sizes
samp_sum<-sample_sums(physeq)
samp_sum_moist<-sample_sums(subset_samples(physeq,soil_moisture =="1" | soil_moisture =="wet"))
samp_sum_dry<-sample_sums(subset_samples(physeq,soil_moisture =="0" | soil_moisture =="dry"))

# generate histograms showing % of contaminating reads
hist(sample_sums(subset_taxa(physeq, Order=="o__Chloroplast"))/samp_sum, breaks = 100)
hist(sample_sums(subset_taxa(physeq, Family=="f__Mitochondria"))/samp_sum, breaks = 100)

hist(sample_sums(subset_taxa(subset_samples(physeq,soil_moisture =="1" | soil_moisture =="wet"), Order=="o__Chloroplast"))/samp_sum_moist, breaks = 100)
hist(sample_sums(subset_taxa(subset_samples(physeq,soil_moisture =="1" | soil_moisture =="wet"), Family=="f__Mitochondria"))/samp_sum_moist, breaks = 100)

hist(sample_sums(subset_taxa(subset_samples(physeq,soil_moisture =="0" | soil_moisture =="dry"), Order=="o__Chloroplast"))/samp_sum_dry, breaks = 100)
hist(sample_sums(subset_taxa(subset_samples(physeq,soil_moisture =="0" | soil_moisture =="dry"), Family=="f__Mitochondria"))/samp_sum_dry, breaks = 100)







#load and run  the fucntion that will remove the plant DNA. press F2 after selecting the custom function to open it on a new tab
physeq<- remove_Chloroplast_Mitochondria(physeq)

# This will check if you still have those taxa in your phyloseq object. if the output is FALSE, then you got rid of them
"o__Chloroplast" %in% tax_table(physeq_clean)
"f__Mitochondria" %in% tax_table(physeq_clean)













######## ----- ITS data ---------- ########

##################### load taxonomy based on full unite reference
physeq_all<-import_biom("./data/Amplicon_sequencing/ITS_phyloseq_input_FeatureTable_metadata_taxonomy.biom",
                        refseqfilename="./data/Amplicon_sequencing/ITS_dna-sequences.fasta")

# removes the horrible ASV naming from qiime2 into ASV numbers
physeq_all<-backup_and_rename(physeq_all) 
########################

#Check how many are fungi and plant ITS

ntaxa(subset_taxa(raw_fun_ps, Kingdom=="k__Fungi"))/ntaxa(raw_fun_ps)*100
ntaxa(subset_taxa(raw_fun_ps, Kingdom=="k__Viridiplantae"))/ntaxa(raw_fun_ps)*100
ntaxa(subset_taxa(raw_fun_ps, Kingdom=="k__Rhizaria"))/ntaxa(raw_fun_ps)*100
ntaxa(subset_taxa(raw_fun_ps, is.na(raw_fun_ps@tax_table)[,1]))/ntaxa(raw_fun_ps)*100


sum(sample_sums(subset_taxa(raw_fun_ps, Domain=="k__Fungi")))/sum(sample_sums((raw_fun_ps)))


# this checks if we have NA in the identification (not found in reference)
check_ratio_of_non_NA<-function(ps_object){
  
  Domain<-ntaxa(subset_taxa(ps_object, Domain!="NA" & Domain!= "Unassigned"))/ntaxa(ps_object)*100
  Phylum<-ntaxa(subset_taxa(ps_object, Phylum!="NA"))/ntaxa(ps_object)*100
  Class<-ntaxa(subset_taxa(ps_object, Class!="NA"))/ntaxa(ps_object)*100
  Order<-ntaxa(subset_taxa(ps_object, Order!="NA"))/ntaxa(ps_object)*100
  Family<-ntaxa(subset_taxa(ps_object, Family!="NA"))/ntaxa(ps_object)*100
  Genus<-ntaxa(subset_taxa(ps_object, Genus!="NA"))/ntaxa(ps_object)*100
  Species<-ntaxa(subset_taxa(ps_object, Species!="NA"))/ntaxa(ps_object)*100
  
  output<-list(Domain,Phylum,Class,Order,Family,Genus,Species)
  names(output)<-c("Domain","Phylum","Class","Order","Family","Genus","Species")
  return(output)
}

check_ratio_of_non_NA(raw_fun_ps)

# this checks identified taxa (not NA  or unidentified)
check_ratio_of_identified<-function(ps_object){
  Domain<-ntaxa(subset_taxa(ps_object, Domain!="NA"& Domain!="k__unidentified" & Domain!= "Unassigned"))/ntaxa(ps_object)*100
  Phylum<-ntaxa(subset_taxa(ps_object, Phylum!="NA"& Phylum!="p__unidentified"))/ntaxa(ps_object)*100
  Class<-ntaxa(subset_taxa(ps_object, Class!="NA" & Class!="c__unidentified"))/ntaxa(ps_object)*100
  Order<-ntaxa(subset_taxa(ps_object, Order!="NA"& Order!="o__unidentified"))/ntaxa(ps_object)*100
  Family<-ntaxa(subset_taxa(ps_object, Family!="NA"& Family!="f__unidentified"))/ntaxa(ps_object)*100
  Genus<-ntaxa(subset_taxa(ps_object, Genus!="NA"& Genus!="g__unidentified"))/ntaxa(ps_object)*100
  Species<-ntaxa(subset_taxa(ps_object, Species!="NA"& Species!="s__unidentified"))/ntaxa(ps_object)*100
  
  output<-list(Domain,Phylum,Class,Order,Family,Genus,Species)
  names(output)<-c("Domain","Phylum","Class","Order","Family","Genus","Species")
  return(output)
}


check_ratio_of_identified(raw_fun_ps)
check_ratio_of_identified(subset_taxa(raw_fun_ps, Domain=="k__Fungi"))