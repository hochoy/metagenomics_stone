# Analyze diamond alignments

## 0. set working directory, set path, make output directory
setwd("~/Analysis_of_gut/metagenomics/")
Sys.setenv(PATH = readLines("scripts_metagenomics/0_setup/bash_path.txt"))
system("mkdir ~/Analysis_of_gut/metagenomics/output/6_blast_for_oxalate_oxidoreductase/")


## 1. load libraries
library(R.utils)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
source("./scripts_metagenomics/0_setup/0_functions.R")


## 2. read in mapping file
oxoxr_db_header <- read_lines("./output/6_blast_for_oxalate_oxidoreductase/diamond_db/oxalateoxidoreductase_header.fasta") 

oxoxr_db_map <- data.frame(header = oxoxr_db_header) %>% 
  mutate(header = as.character(header),
         prot_id = str_extract(header,">[^ ]+ ") %>% trimws(),
         description = str_replace(header,">[^ ]+ ","") %>% trimws(),
         nchar_id = nchar(prot_id)) %>% 
  mutate(subject = prot_id %>% str_replace(">",""))

## 3. read in read depth
readdepth_df <- read.csv("./output/2_qc_and_seqdepth/readdepth_df.csv",header = T)

## 4. identify files in the alignment output folder
aligned_files_df <- data.frame(aligned_files = list.files("./output/6_blast_for_oxalate_oxidoreductase/diamond_alignment/",
                                                          pattern = "*alignment*",
                                                          full.names = T)) %>% 
  mutate(sampleName = str_extract(aligned_files,"CTRLA[0-9AB]{1,4}|PT[0-9AB]{1,4}") %>% trimws()) %>% 
  mutate(aligned_files = as.character(aligned_files))

## 5. function to read in diamond blastx output
get_genecount_from_alignment <- function(aligned_filepath,sampleName) {
  df_raw <- read_delim(aligned_filepath,
             col_names = F,
             delim="\t",
             quote = "",
             comment = "") %>% 
    mutate(sampleName = sampleName) %>% 
    select(sampleName,everything())
  
  new_colnames <- c("sampleName","query","subject",
                    "perc_identity","aln_len","n_mismatch","n_gaps",
                    "query_aln_start","query_aln_end","subj_aln_start","subj_aln_end",
                    "evalue","bitscore")
  df_renamed <- df_raw %>% `colnames<-`(new_colnames)
  
  df_described <- df_renamed %>% 
    inner_join(oxoxr_db_map %>% 
                 select(subject,description),by="subject") %>% 
    mutate(gene_name = str_extract(description,"Oxalate oxidoreductase subunit (alpha|beta|delta)") %>% trimws()) %>% 
    select(sampleName,query,subject,description,gene_name,everything())
  
  # initial diamond alignment filtered for evalue < 0.001; no filter for bitscore; no filter for identity
  df_best_hit <- df_described %>% 
    rbind(df_described[1,]) %>% 
    group_by(query) %>% 
    top_n(1,bitscore) %>%  # filter for best bitscore hit for each query
    slice(1) %>% # if top_n has any draws, keep only 1 hit
    ungroup()
  
  df_genecount <- df_best_hit %>% 
    group_by(sampleName,gene_name) %>% 
    summarize(count = n())
  
  df_genecount_normalized <- df_genecount %>% 
    left_join(readdepth_df,by="sampleName") %>% 
    mutate(count_normalized = count/n_reads * 1000000)
  
  list(df_genecount_normalized)
}

# extract genecount from alignment for each sample
oxr_genecount <- mapply(get_genecount_from_alignment,
                                aligned_filepath = aligned_files_df$aligned_files,
                                sampleName = aligned_files_df$sampleName) %>% 
  do.call(what=rbind)

write_rds(oxr_genecount,
          "./output/6_blast_for_oxalate_oxidoreductase/oxr_genecount.rds")

# keep only normalized count for analysis. reshape for analysis
oxr_normalized <- oxr_genecount %>% 
  select(sampleName,gene_name,count_normalized) %>% 
  spread(gene_name,count_normalized) 

oxr_normalized[is.na(oxr_normalized)] <- 0
  
oxr_normalized <- oxr_normalized %>% 
  gather("gene_name","count_normalized",matches("Oxalate oxidoreductase"))

# read in metadata  
metadata_final <- read.csv("./output/1_metadata/metadata_final.csv",header = T)
metadata_basic <- metadata_final %>% select(sampleName,stoneFormer,pair)

# join oxidoreductase data to metadata
oxr_meta <- oxr_normalized %>% 
  left_join(metadata_basic,by="sampleName")

## calculate oxr presence
oxr_presence <- oxr_meta %>% 
  mutate(presence = ifelse(count_normalized > 0,1,0)) %>% 
  group_by(gene_name,stoneFormer) %>% 
  summarize(presence = sum(presence)) %>% 
  spread(stoneFormer,presence)

# compare abundance of oxidoreductase between pt and ctrl
oxr_abundance <- oxr_meta %>% 
  arrange(pair) %>% 
  compare_counts(count_column = "count_normalized",
                 grouping_column = "gene_name") %>% 
  filter(test == "wilcoxon_signed_ranked") %>% 
  mutate(pval = round(pval,3)) %>% 
  mutate(mean_control = sapply(control,function(x){x %>% str_split(",") %>% unlist() %>% as.numeric() %>% mean() %>% round(1)}),
         sd_control = sapply(control,function(x){x %>% str_split(",") %>% unlist() %>% as.numeric() %>% sd() %>% round(1)}),
         mean_patient = sapply(patient,function(x){x %>% str_split(",") %>% unlist() %>% as.numeric() %>% mean() %>% round(1)}),
         sd_patient = sapply(patient,function(x){x %>% str_split(",") %>% unlist() %>% as.numeric() %>% sd() %>% round(1)}))


# join abundance comparison to presence and write out the table
oxr_presence_abundance <- oxr_presence %>% 
  left_join(oxr_abundance %>% select(gene_name,pval,higher,mean_patient,sd_patient,mean_control,sd_control),by="gene_name") %>% 
  select(gene_name,
         control,patient,
         mean_patient,sd_patient,mean_control,sd_control,
         pval,
         higher)

xlsx::write.xlsx(oxr_presence_abundance %>% as.data.frame(),
                 "./output/5_oxalate_degrading_genes/oxr_presence_abundance.xlsx",col.names = T,row.names = F)

# compare relative abundance of subunits of oxalate oxidoreductase with each other
oxr_forplot <- oxr_genecount %>% 
  left_join(metadata_basic,by="sampleName") %>% 
  gather("count_type","count",matches("count|count_normalized"))  %>% 
  group_by(count_type,pair,gene_name) %>% 
  
  mutate(compare_count = paste0(count,collapse=","),
         compare_grp = paste0(stoneFormer,collapse=","),
         compare_pair = paste0(pair,collapse=",")) %>% 
  mutate(higher = mapply(which_higher,
                         comma_separated_num = compare_count,
                         comma_separated_grp = compare_grp)) %>% 
  mutate(gene_descr = str_wrap(gene_name,width = 20)) %>% 
  ungroup()

# (oxr_abscount_boxplot <- oxr_forplot %>% 
#     filter(count_type == "count") %>% 
#     ggplot(aes(x=stoneFormer,y=count,fill=stoneFormer)) +
#     geom_boxplot() +
#     facet_wrap(~gene_descr,scales="free_y") +
#     scale_fill_manual(values=c("#00BFC4","#F8766D"))
# )

(oxr_normalized_count_boxplot <- oxr_forplot %>% 
    filter(count_type == "count_normalized") %>% 
    ggplot(aes(x=stoneFormer,y=count,fill=stoneFormer)) +
    geom_boxplot() +
    facet_wrap(~gene_descr,scales="free_y") +
    scale_fill_manual(values=c("#00BFC4","#F8766D")) +
    ylab("Gene count per million reads") +
    xlab("")
)

## plot for thesis
ggsave(filename = "./output/6_blast_for_oxalate_oxidoreductase/oxr_boxplot.png",plot = oxr_normalized_count_boxplot,
       width = 25,height=10,units = "cm",dpi = 350)

(oxr_abscount_pairplot <- oxr_forplot %>% 
  filter(count_type == "count") %>% 
  ggplot(aes(x=stoneFormer,y=count)) +
  geom_point(aes(color=stoneFormer)) +
  facet_wrap(~ gene_descr) +
  geom_line(aes(group = pair,color=higher)) +
  scale_color_manual(values=c("#00BFC4","black","#F8766D")))
  
(oxr_normcount_pairplot <- oxr_forplot %>% 
    filter(count_type == "count_normalized") %>% 
    ggplot(aes(x=stoneFormer,y=count)) +
    geom_point(aes(color=stoneFormer)) +
    facet_wrap(~ gene_descr) +
    geom_line(aes(group = pair,color=higher)) +
    scale_color_manual(values=c("#00BFC4","black","#F8766D")) +
    ylab("gene count per million reads"))


# plot out image for results


