# oxalate degrading genes

## 0. set working directory, set path, make output directory
setwd("~/Analysis_of_gut/metagenomics/")
Sys.setenv(PATH = readLines("script_metagenomics/0_setup/bash_path.txt"))
system("mkdir ~/Analysis_of_gut/metagenomics/output/5_oxalate_degrading_genes")

## 1. load required libraries
library(coin)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(readr)
source("./script_metagenomics/0_setup/0_functions.R")


## 2. read in the kegg tpm results
kegg_tpm_results <- read_tsv("./input/tpm_results/kegg_tpm_results.txt",col_names = T)
kegg_tpm_results[is.na(kegg_tpm_results)] <- 0


## 3. read in the kegg hierarchy + extract enzyme or transporter commission id
## result: there are 8077 unique kegg orthology ids in our kegg hierarchy table
kegg_hierarchy <- read_tsv("./input/kegg/kegg_hierarchy.txt")
kegg_hierarchy_ec <- kegg_hierarchy %>% 
  mutate(EC = str_match(gene_descr,"[ET]{1}C:([a-z,A-Z,0-9-]+\\.[a-z,A-Z,0-9-]+\\.[a-z,A-Z,0-9-]+\\.[a-z,A-Z,0-9-]+)")[,2] %>% trimws())
kegg_hierarchy_ec$KO %>% unique() %>% length() # 8077 


## 4. remove duplicate KEGG ids
## there are duplicates because each KEGG orthology id might belong to multiple metabolic pathways
kegg_hierarchy_dedup <- kegg_hierarchy_ec %>% 
  select(-hierarchy1, # remove the upper hierarchies to deduplicate kegg id
         -hierarchy2,
         -pathway,
         -pathway_name) %>%  
  unique()


## 5. identify which of the 8077 kegg ids are present in our tpm results
## result: 2816 of the 8077 KEGG orthology ids (34.8%) were detected in at least 1 of our samples
tpm_hierarchized <- kegg_tpm_results %>% 
  inner_join(kegg_hierarchy_dedup,
             by="KO") %>% 
  select(KO, 
         gene_descr,
         EC, 
         everything())
tpm_hierarchized$KO %>% unique() %>% length() # 2816


## 6. compile oxalate-associated enzymes
## Done by searching the KEGG website for genes that have oxalate, formate, formyl-coa, oxalyl-coa as substrate, product or both
kegg_oxalate_enzymes <- read.csv("./input/kegg/kegg_oxalate_enzymes.csv",header = T,na.strings = c("","NA"))
bacterial_oxalate_enzymes <- kegg_oxalate_enzymes %>% filter(species == "bacterial") %>% 
  mutate(enzyme_index = 1:nrow(.))

## 7. filter our data for only oxalate enzymes
## by joining our data to the oxalate enzymes via kegg id
oxalate_enzyme_tpm <- bacterial_oxalate_enzymes %>% 
  select(-species,-other_names,-source,-kegg_link,-other_link) %>% 
  left_join(tpm_hierarchized %>% select(-EC),
            by="KO")

## 8. add in metadata
metadata_final <- read.csv("./output/1_metadata/metadata_final.csv",header = T)
metadata_basic <- metadata_final %>% select(sampleName,stoneFormer,pair)

oxalate_enzyme_meta <- oxalate_enzyme_tpm %>% 
  gather("sampleName","tpm",matches("CTRLA[0-9,A-Z]+|PT[0-9,A-Z]+")) %>% 
  inner_join(metadata_basic, by="sampleName")

## 9. calculate presence
oxalate_enzyme_presence <- oxalate_enzyme_meta %>% 
  mutate(presence = ifelse(tpm > 0,1,0)) %>% 
  select(-sampleName,-tpm,-pair) %>% 
  group_by_at(vars(-presence)) %>% 
  summarize(presence = sum(presence)) %>% 
  spread(stoneFormer,presence) %>% 
  rename(n_patient = patient,
         n_control = control) %>% 
  ungroup()

## 10. compare abundance
oxalate_enzyme_abundance <- oxalate_enzyme_meta %>% 
  filter(!is.na(tpm)) %>% # remove non-detected enzymes
  select(enzyme_index,sampleName,tpm,stoneFormer,pair) %>% 
  compare_counts(count_column = "tpm",
                 grouping_column = "enzyme_index") %>% 
  filter(test == "wilcoxon_signed_ranked") %>% 
  mutate(patient = str_split(patient,",") %>% lapply(as.numeric),
         control = str_split(control,",") %>% lapply(as.numeric))

## 11. add in abundance comparison into presence df
oxalate_enzyme_abundance_presence <- oxalate_enzyme_abundance %>% 
  left_join(oxalate_enzyme_presence %>% select(enzyme_index,n_patient,n_control),by="enzyme_index")

### 10b.1 Add all oxalate enzymes together to compare
oxalate_enzyme_sum <- oxalate_enzyme_meta %>% 
  group_by(stoneFormer,pair,sampleName) %>% 
  summarize(sum_tpm = sum(tpm,na.rm = T)) %>% 
  ungroup() %>% 
  mutate(stoneFormer = factor(stoneFormer,levels=c("patient","control"))) 

### 10b.2 compare
oxalate_enzyme_sum_test <- oxalate_enzyme_sum %>% 
  mutate(enzyme_index = 999) %>% 
  compare_counts("sum_tpm",grouping_column="enzyme_index") %>% 
  mutate(control= str_split(control,",") %>% lapply(as.numeric),
         patient= str_split(patient,",") %>% lapply(as.numeric),
         n_control = lapply(control, function(vec){length(which(as.numeric(vec) != 0))}),
         n_patient = lapply(patient, function(vec){length(which(as.numeric(vec) != 0))})
  ) %>% 
  filter(test == "wilcoxon_signed_ranked")


### 10b.3 combine individual abundance with sum abundance
oxalate_enzyme_alltests <- rbind(oxalate_enzyme_abundance_presence,
                                 oxalate_enzyme_sum_test)

### 11. extract mean and se
oxalate_enzyme_fortable <- oxalate_enzyme_alltests %>% 
  mutate(patient_mean = sapply(patient,function(x){mean(x) %>% round(1)}),
         control_mean = sapply(control,function(x){mean(x) %>% round(1)}),
         patient_se = sapply(patient,function(x){se(x) %>% round(1)}),
         control_se = sapply(control,function(x){se(x) %>% round(1)}),
         patient_label = paste(patient_mean,plusminus,patient_se),
         control_label = paste(control_mean,plusminus,control_se))  

## 12. Add back enzyme name
oxalate_enzyme_genename <- oxalate_enzyme_fortable %>% 
  full_join(bacterial_oxalate_enzymes %>% select(gene_name,enzyme_index,gene_abbrev,KO,EC),by="enzyme_index") %>% 
  mutate(gene_name = as.character(gene_name),
         KO=as.character(KO),
         EC=as.character(EC)) %>% 
  replace_na(list(gene_name="All oxalate-associated enzymes"))

## 12. Draw table
oxalate_enzyme_table <- oxalate_enzyme_genename %>% 
  arrange(enzyme_index) %>% 
  select(gene_name,
         n_patient,n_control,
         patient_label,control_label,
         pval,higher,
         gene_abbrev,
         KO,EC
  ) %>% 
  mutate(gene_abbrev = as.character(gene_abbrev)) %>% 
  mutate(pval = round(pval,3)) %>% 
  rename(`Gene name`=gene_name,
         `Patient (n)`=n_patient,`Control (n)`=n_control,
         `Patient (tpm)`=patient_label,
         `Control (tpm)`=control_label,
         `Gene abbreviation`=gene_abbrev) 

oxalate_enzyme_table[is.na(oxalate_enzyme_table)] <- "-"
oxalate_enzyme_table$`Patient (n)`[sapply(FUN=is.null,oxalate_enzyme_table$`Patient (n)`)] <- "-"
oxalate_enzyme_table$`Control (n)`[sapply(FUN=is.null,oxalate_enzyme_table$`Control (n)`)] <- "-"

xlsx::write.xlsx(oxalate_enzyme_table %>% as.data.frame(),
                 "./output/5_oxalate_degrading_genes/oxalate_enzyme_table.xlsx",col.names = T,row.names = F)
# table for thesis
xlsx::write.xlsx(oxalate_enzyme_table %>% as.data.frame(),
                 "../thesis/oxalate_enzyme_c2.xlsx",col.names = T,row.names = F)


## 13. Draw separate enzyme reaction table
oxalate_reaction_table <- bacterial_oxalate_enzymes %>% 
  select(gene_name,
         reaction,
         gene_abbrev,
         KO, EC
  )%>% 
  mutate(gene_abbrev = as.character(gene_abbrev),
         KO = as.character(KO))
oxalate_reaction_table[is.na(oxalate_reaction_table)] <- "-"

xlsx::write.xlsx(oxalate_reaction_table %>% as.data.frame(),
                 "./output/5_oxalate_degrading_genes/oxalate_reaction_table.xlsx",col.names = T,row.names = F)

# table for thesis
xlsx::write.xlsx(oxalate_reaction_table %>% as.data.frame(),
                 "../thesis/oxalate_reaction_c2.xlsx",col.names = T,row.names = F)



