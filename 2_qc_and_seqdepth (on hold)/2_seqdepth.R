# Analyze sequencing depth, orf, and KO

## 0. set working directory, set path, make output directory
setwd("~/Analysis_of_gut/metagenomics/")
Sys.setenv(PATH = readLines("scripts_metagenomics/0_setup/bash_path.txt"))
system("mkdir ~/Analysis_of_gut/metagenomics/output/2_qc_and_seqdepth/")

## 1. load library
library(dplyr)
library(stringr)
source("./scripts_metagenomics/0_setup/0_functions.R")

# FOR THESIS

# Read in reads before and after kneading (trimmomatic + contaminant removal)
rawread_firstseq_checkknead <- read_rds("~/aNALYSIS_PIPELINE/data/evaluation/ME1_prep_rawreads_humann2/rawread_firstseq_checkknead.rds")
rawread_secondseq_checkknead <- read_rds("~/aNALYSIS_PIPELINE/data/evaluation/ME1_prep_rawreads_humann2/rawread_secondseq_checkknead.rds")
rawread_thirdseq_checkknead <- read_rds("~/aNALYSIS_PIPELINE/data/evaluation/ME1_prep_rawreads_humann2/rawread_thirdseq_checkknead.rds")

# Calculate original reads and final reads
rawread_all_checkknead <- rbind(rawread_firstseq_checkknead,
                                rawread_secondseq_checkknead,
                                rawread_thirdseq_checkknead)

before_after_knead <- rawread_all_checkknead %>% 
  mutate(sum_kneaded = sum(knead_paired1_reads,
                           knead_paired2_reads,
                           knead_unmatched1_reads,
                           knead_unmatched2_reads),
         total_reads = sum(forward_reads,
                           reverse_reads)) %>% 
  mutate(sampleGroup = ifelse(str_detect(sampleName,"CTRL"),"CTRL","PT"))

# add metadata
metadata_final <- read.csv("~/Analysis_of_gut/metagenomics/output/1_metadata/metadata_final.csv")

basic_beforeafterknead <- before_after_knead %>% 
  left_join(metadata_final %>% select(sampleName,stoneFormer,pair,XXX_sequencing_round)) %>% 
  select(sampleName,stoneFormer,pair,XXX_sequencing_round,total_reads,sum_kneaded) %>% 
  rename(sampleName = sampleName,
         `Sample group` = stoneFormer,
         `Pair` = pair,
         `Sequencing batch` = XXX_sequencing_round,
         `Total sequencing reads` = total_reads,
         `Trimmed reads` = sum_kneaded) %>% 
  arrange(Pair,`Sample group`)

# add orf and ko results
orf_stats_results <- read.delim("~/aNALYSIS_PIPELINE/data/metagenomics/MC3_prep_rpkmoutput/orf_stats_results.txt",
                                header = T,sep = "\t") %>% 
  mutate(sampleName = toupper(sampleName))

ORF_df <- orf_stats_results %>% 
  group_by(sampleName) %>% 
  summarize(total_orfs = sum(n_orfs),
            total_tpm = sum(sum_tpm))

KO_valid_df <- orf_stats_results %>% 
  filter(str_detect(status,"has [12]{1}")) %>% 
  group_by(sampleName) %>% 
  summarize(total_orfs_with_KO = sum(n_orfs),
            total_tpm_with_KO = sum(sum_tpm))

orf_ko_df <- ORF_df %>% 
  left_join(KO_valid_df) %>% 
  select(sampleName,
         total_orfs,
         total_orfs_with_KO,
         total_tpm,
         total_tpm_with_KO)

ko_tpm_table <- basic_beforeafterknead %>% 
  left_join(orf_ko_df) %>% 
  rename(`Sample name` = sampleName,
         `Total ORFs (predicted by Prodigal)` = total_orfs,
         `Total ORFs with valid KEGG Orthology IDs` = total_orfs_with_KO,
         `Total TPM` = total_tpm,
         `Total TPM for ORFs with valid KEGG Orthology IDs` = total_tpm_with_KO)

ko_tpm_table_writeout <- ko_tpm_table %>% 
  mutate(`Total sequencing reads` = formatC(`Total sequencing reads`, format="d", big.mark=","),
         `Trimmed reads` = formatC(`Trimmed reads`, format="d", big.mark=","),
         `Total ORFs (predicted by Prodigal)` = formatC(`Total ORFs (predicted by Prodigal)`, format="d", big.mark=","),
         `Total ORFs with valid KEGG Orthology IDs` = formatC(`Total ORFs with valid KEGG Orthology IDs`, format="d", big.mark=","),
         `Total TPM for ORFs with valid KEGG Orthology IDs` = formatC(`Total TPM for ORFs with valid KEGG Orthology IDs`, format="d", big.mark=","))
write.csv(ko_tpm_table_writeout,"./output/2_qc_and_seqdepth/ko_tpm_table.csv",row.names = F)

# draw contrasted plot
ko_tpm_table %>% 
  group_by(Pair) %>% 
  mutate(ratio_sequencingreads = max(`Total sequencing reads`)/min(`Total sequencing reads`),
         ratio_trimmedreads = max(`Trimmed reads`)/min(`Trimmed reads`),
         ratio_orfs = max(`Total ORFs (predicted by Prodigal)`)/min(`Total ORFs (predicted by Prodigal)`),
         ratio_valid_orfs = max(`Total ORFs with valid KEGG Orthology IDs`)/min(`Total ORFs with valid KEGG Orthology IDs`)) %>% View()

annot_font_size = 4

(metagenomics_readdepth_ctplot <- ko_tpm_table %>% 
    mutate(color = ifelse(Pair != 1,"grey",ifelse(`Sample group` == "patient","red","blue"))) %>%
    mutate(color = factor(color,levels = c("blue","red","grey"))) %>% 
    ggplot(aes(x=Pair,y=`Total sequencing reads`,group=`Sample group`,fill=color)) +
    geom_bar(stat="identity",position = position_dodge(width=0.8),width = 0.7) +
    
    scale_fill_manual(labels = c("Control","Patient","<3-fold difference"),values=c("blue" = "#00BFC4","red" = "#F8766D","grey" = "grey")) +
    scale_y_continuous(label=scales::comma,breaks = seq(0,180000000,50000000))+
    scale_x_continuous(breaks=seq(1,17,1)) +
    xlab("sample pair") +
    ylab("total sequencing reads") +
    labs(fill="group")+
    # annotate("text",x=7,y=30000,label = "12.8 X",color = "#00BFC4",size =annot_font_size) +
    annotate("text",x=1,y=155000000,label = "3.7 X",color = "#F8766D",size =annot_font_size) +
    theme1
)

ggsave(filename = "./output/2_qc_and_seqdepth/metagenomics_readdepth_ctplot.tiff",
       plot = metagenomics_readdepth_ctplot,
       dpi = 280,
       width = 20,
       height = 14,
       units = "cm")

(metagenomics_trimmed_readdepth_ctplot <- ko_tpm_table %>% 
    mutate(color = ifelse(Pair != 1,"grey",ifelse(`Sample group` == "patient","red","blue"))) %>%
    mutate(color = factor(color,levels = c("blue","red","grey"))) %>% 
    ggplot(aes(x=Pair,y=`Trimmed reads`,group=`Sample group`,fill=color)) +
    geom_bar(stat="identity",position = position_dodge(width=0.8),width = 0.7) +
    
    scale_fill_manual(labels = c("Control","Patient","<3-fold difference"),values=c("blue" = "#00BFC4","red" = "#F8766D","grey" = "grey")) +
    scale_y_continuous(label=scales::comma,breaks = seq(0,140000000,25000000))+
    scale_x_continuous(breaks=seq(1,17,1)) +
    xlab("sample pair") +
    ylab("total trimmed sequencing reads") +
    labs(fill="group")+
    # annotate("text",x=7,y=30000,label = "12.8 X",color = "#00BFC4",size =annot_font_size) +
    annotate("text",x=1,y=130000000,label = "3.8 X",color = "#F8766D",size =annot_font_size) +
    theme1
)

ggsave(filename = "./output/2_qc_and_seqdepth/metagenomics_trimmed_readdepth_ctplot.tiff",
       plot = metagenomics_trimmed_readdepth_ctplot,
       dpi = 280,
       width = 20,
       height = 14,
       units = "cm")

(orf_depth_ctplot <- ko_tpm_table %>% 
    mutate(color = "grey") %>% 
    ggplot(aes(x=Pair,y=`Total ORFs (predicted by Prodigal)`,group=`Sample group`,fill=color)) +
    geom_bar(stat="identity",position = position_dodge(width=0.8),width = 0.7) +
    
    scale_fill_manual(labels = c("<3-fold difference"),values=c("blue" = "#00BFC4","red" = "#F8766D","grey" = "grey")) +
    scale_y_continuous(label=scales::comma,breaks = seq(0,600000,100000))+
    scale_x_continuous(breaks=seq(1,17,1)) +
    xlab("sample pair") +
    ylab("total ORFs (predicted by Prodigal)") +
    labs(fill="group")+
    # annotate("text",x=7,y=30000,label = "12.8 X",color = "#00BFC4",size =annot_font_size) +
    # annotate("text",x=1,y=155000000,label = "3.7 X",color = "#F8766D",size =annot_font_size) +
    theme1
)

ggsave(filename = "./output/2_qc_and_seqdepth/orf_depth_ctplot.tiff",
       plot = orf_depth_ctplot,
       dpi = 280,
       width = 20,
       height = 14,
       units = "cm")


(orf_with_ko_depth_ctplot <- ko_tpm_table %>% 
    mutate(color = "grey") %>% 
    ggplot(aes(x=Pair,y=`Total ORFs with valid KEGG Orthology IDs`,group=`Sample group`,fill=color)) +
    geom_bar(stat="identity",position = position_dodge(width=0.8),width = 0.7) +
    
    scale_fill_manual(labels = c("<3-fold difference"),values=c("blue" = "#00BFC4","red" = "#F8766D","grey" = "grey")) +
    scale_y_continuous(label=scales::comma,breaks = seq(0,600000,100000))+
    scale_x_continuous(breaks=seq(1,17,1)) +
    xlab("sample pair") +
    ylab("Total ORFs with valid KEGG Orthology IDs") +
    labs(fill="group")+
    # annotate("text",x=7,y=30000,label = "12.8 X",color = "#00BFC4",size =annot_font_size) +
    # annotate("text",x=1,y=155000000,label = "3.7 X",color = "#F8766D",size =annot_font_size) +
    theme1
)

ggsave(filename = "./output/2_qc_and_seqdepth/orf_with_ko_depth_ctplot.tiff",
       plot = orf_with_ko_depth_ctplot,
       dpi = 280,
       width = 20,
       height = 14,
       units = "cm")

ko_tpm_table %>% 
  ggplot(aes(x=`Total sequencing reads`,y=`Total ORFs with valid KEGG Orthology IDs`)) +
  geom_point() +
  geom_smooth(method = "lm")

ko_tpm_table %>% 
  ggplot(aes(x=`Total sequencing reads`,y=`Total ORFs (predicted by Prodigal)`)) +
  geom_point() +
  geom_smooth(method = "lm")

# OLD
# ## 2. find and organize the clean reads
# clean_reads_df <- data.frame(cleanread_files = list.files("./input/clean_reads/",
#                                                           pattern = "merged.gz",
#                                                           full.names = T)) %>% 
#   mutate(sampleName = str_extract(cleanread_files,"CTRLA[0-9AB]{1,4}|PT[0-9AB]{1,4}") %>% trimws())
# 
# ## 3. count the reads in each sample's read file
# readdepth_df <- clean_reads_df %>% 
#   mutate(n_reads = sapply(cleanread_files,
#                           count_gz_seqs)) 
# 
# ## 4. save
# # write.csv(readdepth_df,"./output/2_qc_and_seqdepth/readdepth_df.csv",
# #           row.names = F,col.names = T)
# readdepth_df <- read.csv("./output/2_qc_and_seqdepth/readdepth_df.csv") %>% 
#   mutate(single_end_reads = n_reads/2)