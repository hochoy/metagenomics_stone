# uric acid degrading genes

## 0. set working directory, set path, make output directory
setwd("~/Analysis_of_gut/metagenomics/")
Sys.setenv(PATH = readLines("scripts_metagenomics/0_setup/bash_path.txt"))
system("mkdir ~/Analysis_of_gut/metagenomics/output/5_oxalate_degrading_genes")

## 1. load required libraries
library(coin)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(readr)
source("./scripts_metagenomics/0_setup/0_functions.R")


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


## 6. compile uric-associated enzymes
# ADD:
# Uric acid:
# KEGG pathway: http://www.genome.jp/kegg-bin/show_pathway?map00230+C00366
# KEGG enzymes: http://www.genome.jp/dbget-bin/get_linkdb?-t+enzyme+cpd:C00366 
# enzyme path: https://www.researchgate.net/figure/While-degradation-of-purines-to-uric-acid-is-generally-conserved-among-organisms-the-end_fig1_236693470


## Done by searching the KEGG website for genes that are connected to uric acid directly or to its degradation pathway
kegg_uric_enzymes <- read.csv("./input/kegg/kegg_uric_enzymes.csv",header = T,na.strings = c("","NA"),stringsAsFactors = F)
bacterial_uric_enzymes <- kegg_uric_enzymes %>%
  mutate(enzyme_index = 1:nrow(.))

bacterial_uric_enzymes$other_names[is.na(bacterial_uric_enzymes$other_names)] <- ""

## 7. filter our data for only uric enzymes
## by joining our data to the uric enzymes via kegg id
uric_enzyme_tpm <- bacterial_uric_enzymes %>% 
  select(-species,-other_names,-source,-kegg_link,-other_link) %>% 
  left_join(tpm_hierarchized %>% select(-EC),
            by="KO")

## 8. add in metadata
metadata_final <- read.csv("./output/1_metadata/metadata_final.csv",header = T)
metadata_basic <- metadata_final %>% select(sampleName,stoneFormer,pair)

uric_enzyme_meta <- uric_enzyme_tpm %>% 
  gather("sampleName","tpm",matches("CTRLA[0-9,A-Z]+|PT[0-9,A-Z]+")) %>% 
  inner_join(metadata_basic, by="sampleName")

## 9. calculate presence
uric_enzyme_presence <- uric_enzyme_meta %>% 
  mutate(presence = ifelse(tpm > 0,1,0)) %>% 
  select(-sampleName,-tpm,-pair) %>% 
  group_by_at(vars(-presence)) %>% 
  summarize(presence = sum(presence)) %>% 
  spread(stoneFormer,presence) %>% 
  rename(n_patient = patient,
         n_control = control) %>% 
  ungroup()

## 10. compare abundance
uric_enzyme_abundance <- uric_enzyme_meta %>% 
  filter(!is.na(tpm)) %>% # remove non-detected enzymes
  select(enzyme_index,sampleName,tpm,stoneFormer,pair) %>% 
  compare_counts(count_column = "tpm",
                 grouping_column = "enzyme_index") %>% 
  filter(test == "wilcoxon_signed_ranked") %>% 
  mutate(patient = str_split(patient,",") %>% lapply(as.numeric),
         control = str_split(control,",") %>% lapply(as.numeric))

## 11. add in abundance comparison into presence df
uric_enzyme_abundance_presence <- uric_enzyme_abundance %>% 
  left_join(uric_enzyme_presence %>% select(enzyme_index,n_patient,n_control),by="enzyme_index")

### 10b.1 Add all uric enzymes together to compare
uric_enzyme_sum <- uric_enzyme_meta %>% 
  group_by(stoneFormer,pair,sampleName) %>% 
  summarize(sum_tpm = sum(tpm,na.rm = T)) %>% 
  ungroup() %>% 
  mutate(stoneFormer = factor(stoneFormer,levels=c("patient","control"))) 

### 10b.2 compare
uric_enzyme_sum_test <- uric_enzyme_sum %>% 
  mutate(enzyme_index = 999) %>% 
  compare_counts("sum_tpm",grouping_column="enzyme_index") %>% 
  mutate(control= str_split(control,",") %>% lapply(as.numeric),
         patient= str_split(patient,",") %>% lapply(as.numeric),
         n_control = lapply(control, function(vec){length(which(as.numeric(vec) != 0))}),
         n_patient = lapply(patient, function(vec){length(which(as.numeric(vec) != 0))})
  ) %>% 
  filter(test == "wilcoxon_signed_ranked")


### 10b.3 combine individual abundance with sum abundance
uric_enzyme_alltests <- rbind(uric_enzyme_abundance_presence,
                                 uric_enzyme_sum_test)

### 11. extract mean and se
uric_enzyme_fortable <- uric_enzyme_alltests %>% 
  mutate(patient_mean = sapply(patient,function(x){mean(x) %>% round(1)}),
         control_mean = sapply(control,function(x){mean(x) %>% round(1)}),
         patient_se = sapply(patient,function(x){se(x) %>% round(1)}),
         control_se = sapply(control,function(x){se(x) %>% round(1)}),
         patient_label = paste(patient_mean,plusminus,patient_se),
         control_label = paste(control_mean,plusminus,control_se))  

## 12. Add back enzyme name
uric_enzyme_genename <- uric_enzyme_fortable %>% 
  full_join(bacterial_uric_enzymes %>% select(gene_name,other_names,enzyme_index,gene_abbrev,KO,EC),by="enzyme_index") %>% 
  mutate(gene_name = paste0(gene_name," (",other_names,")",sep=""),
         KO=as.character(KO),
         EC=as.character(EC)) %>% 
  mutate(gene_name = str_replace(gene_name,"NA \\(NA\\)","All uric-acid-degrading enzymes"))

## 12. Draw table
uric_enzyme_table <- uric_enzyme_genename %>% 
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

uric_enzyme_table[is.na(uric_enzyme_table)] <- "-"
uric_enzyme_table$`Patient (n)`[sapply(FUN=is.null,uric_enzyme_table$`Patient (n)`)] <- "-"
uric_enzyme_table$`Control (n)`[sapply(FUN=is.null,uric_enzyme_table$`Control (n)`)] <- "-"

xlsx::write.xlsx(uric_enzyme_table %>% as.data.frame(),
                 "./output/5_oxalate_degrading_genes/uric_enzyme_table.xlsx",col.names = T,row.names = F)
# table for thesis
xlsx::write.xlsx(uric_enzyme_table %>% as.data.frame(),
                 "../thesis/chapter1/uric_enzyme_table_c1.xlsx",col.names = T,row.names = F)


## 13. Draw separate enzyme reaction table
uric_reaction_table <- bacterial_uric_enzymes %>% 
  select(gene_name,
         reaction,
         gene_abbrev,
         KO, EC
  )%>% 
  mutate(gene_abbrev = as.character(gene_abbrev),
         KO = as.character(KO))
uric_reaction_table[is.na(uric_reaction_table)] <- "-"

xlsx::write.xlsx(uric_reaction_table %>% as.data.frame(),
                 "./output/5_oxalate_degrading_genes/uric_reaction_table.xlsx",col.names = T,row.names = F)

# table for thesis
xlsx::write.xlsx(uric_reaction_table %>% as.data.frame(),
                 "../thesis/chapter1/uric_reaction_table_c1.xlsx",col.names = T,row.names = F)






