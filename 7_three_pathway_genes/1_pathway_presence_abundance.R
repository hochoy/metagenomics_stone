# analyzing three pathway genes

## 0. set working directory, set path, make output directory
setwd("~/Analysis_of_gut/metagenomics/")
Sys.setenv(PATH = readLines("./scripts_metagenomics/0_setup/bash_path.txt"))
system("mkdir ~/Analysis_of_gut/metagenomics/output/7_three_pathway_genes")

## 1. load required libraries
library(coin)
library(ggplot2)
library(stringr)
library(tidyr)
library(readr)
library(dplyr)
source("./scripts_metagenomics/0_setup/0_functions.R")

## 2. read in the kegg tpm results
kegg_tpm_results <- read_tsv("./input/tpm_results/kegg_tpm_results.txt",col_names = T)
kegg_tpm_results[is.na(kegg_tpm_results)] <- 0


## 3. read in the kegg hierarchy + extract enzyme or transporter commission id
## result: there are 8077 unique kegg orthology ids in our kegg hierarchy table
# regex for capturing and removing EC
ec_regex <- "[ET]{1}C:([a-z,A-Z,0-9-]+\\.[a-z,A-Z,0-9-]+\\.[a-z,A-Z,0-9-]+\\.[a-z,A-Z,0-9-]+[ ]{0,1})"
ec_removal_front_regex <- "[ET]{1}C{0,1}:{0,1}([a-z,A-Z,0-9-]+\\.[a-z,A-Z,0-9-]+\\.[a-z,A-Z,0-9-]+\\.[a-z,A-Z,0-9-]+[ ,]{0,1})"
ec_removal_back_regex <- paste0("\\[[ ]{0,1}",ec_regex,"+","\\]")

kegg_hierarchy <- read_tsv("./input/kegg/kegg_hierarchy.txt")
kegg_hierarchy_ec <- kegg_hierarchy %>% 
  mutate(EC = str_match(gene_descr,ec_regex)[,2] %>% trimws()) %>% 
  mutate(gene_abbrev = str_extract(gene_descr,"^.+?;") %>% str_replace(ec_removal_front_regex,"") %>% str_replace(";","") %>%  trimws()) %>% 
  mutate(gene_name = str_replace(gene_descr,"^.+?;","") %>% str_replace(ec_removal_back_regex,"") %>% trimws()) %>% 
  mutate(name_nchar = nchar(gene_name),
         abbrev_nchar = nchar(gene_abbrev))

## 4. join tpm results to kegg hierarchy
tpm_hierarchized <- kegg_hierarchy_ec %>% 
  left_join(kegg_tpm_results,by="KO")



# ADD:
# Uric acid:
# KEGG pathway: http://www.genome.jp/kegg-bin/show_pathway?map00230+C00366
# KEGG enzymes: http://www.genome.jp/dbget-bin/get_linkdb?-t+enzyme+cpd:C00366 
# enzyme path: https://www.researchgate.net/figure/While-degradation-of-purines-to-uric-acid-is-generally-conserved-among-organisms-the-end_fig1_236693470

# Cystine: 
# KEGG pathway: http://www.genome.jp/kegg-bin/show_pathway?map00270+C00491
# KEGG enzymes: http://www.genome.jp/dbget-bin/get_linkdb?-t+enzyme+cpd:C00097
# enzyme path: https://www.researchgate.net/publication/6698410_Functional_Analysis_of_luxS_in_the_Probiotic_Strain_Lactobacillus_rhamnosus_GG_Reveals_a_Central_Metabolic_Role_Important_for_Growth_and_Biofilm_Formation?_sg=Uh9C1emxNDA05iosb7aHoVlyMK8ZqMoMCL7Bzij_s5mlb9RcAQqG2LYw0_P6FAHjsOlxAUng5Q
# 2nd link: http://jb.asm.org/content/194/13/3522/F1.expansion.html 
# cool: https://www.nature.com/articles/s41598-018-19920-y
# cys to sulfur https://www.researchgate.net/figure/Pathways-of-cysteine-degradation-to-H2S-Cysteine-desulfhydrase-is-a-key-enzyme-for_fig3_233888742 
# degrading bacteria: http://jb.asm.org/content/199/16/e00117-17.full 
# old, but full review of all mechanims including transport: https://link.springer.com/chapter/10.1007/7171_2006_060 



## 5. filter and split into the 3 metabolic pathways
buta_pwy <- "ko00650"
ascor_pwy <- "ko00053"
glyox_pwy <- "ko00630"
uric_pwy <- "ko00230"
cyst_pwy <- "ko00270"

buta_tpm <- tpm_hierarchized %>% filter(pathway == buta_pwy)
ascor_tpm <- tpm_hierarchized %>% filter(pathway == ascor_pwy)
glyox_tpm <- tpm_hierarchized %>% filter(pathway == glyox_pwy)
uric_tpm <- tpm_hierarchized %>% filter(pathway == uric_pwy)
cyst_tpm <- tpm_hierarchized %>% filter(pathway == cyst_pwy)

buta_tpm$KO %>% duplicated() %>% which() # pass, no duplicates
ascor_tpm$KO %>% duplicated() %>% which() # pass, no duplicates
glyox_tpm$KO %>% duplicated() %>% which() # pass, no duplicates
uric_tpm$KO %>% duplicated() %>% which() # pass, no duplicates
cyst_tpm$KO %>% duplicated() %>% which() # pass, no duplicates

## 6. add metadata
metadata_final <- read.csv("./output/1_metadata/metadata_final.csv",header = T)
metadata_basic <- metadata_final %>% select(sampleName,stoneFormer,pair)

add_metadata <- function(df) {
  df %>% 
    gather("sampleName","tpm",matches("CTRLA[0-9,A-Z]+|PT[0-9,A-Z]+")) %>% 
    inner_join(metadata_basic, by="sampleName")
}
buta_meta <- buta_tpm %>% add_metadata()
ascor_meta <- ascor_tpm %>% add_metadata()
glyox_meta <- glyox_tpm %>% add_metadata()
uric_meta <- uric_tpm %>% add_metadata()
cyst_meta <- cyst_tpm %>% add_metadata()

write_rds(buta_meta,"./output/7_three_pathway_genes/buta_meta.rds")
write_rds(ascor_meta,"./output/7_three_pathway_genes/ascor_meta.rds")
write_rds(glyox_meta,"./output/7_three_pathway_genes/glyox_meta.rds")
write_rds(uric_meta,"./output/7_three_pathway_genes/uric_meta.rds")
write_rds(cyst_meta,"./output/7_three_pathway_genes/cyst_meta.rds")

## 7. calculate individual gene presence
calculate_presence <- function(df) {
  df %>% 
    mutate(presence = ifelse(tpm > 0,1,0)) %>% 
    select(-sampleName,-tpm,-pair,-hierarchy1,-hierarchy2) %>% 
    group_by_at(vars(-presence)) %>% 
    summarize(presence = sum(presence)) %>% 
    spread(stoneFormer,presence) %>% 
    rename(n_patient = patient,
           n_control = control) %>% 
    ungroup()
}
buta_presence <- buta_meta %>% calculate_presence()
ascor_presence <- ascor_meta %>% calculate_presence()
glyox_presence <- glyox_meta %>% calculate_presence()
uric_presence <- uric_meta %>% calculate_presence()
cyst_presence <- cyst_meta %>% calculate_presence()

## 8. compare individual gene abundance
compare_abundance <- function(df) {
  df %>% 
    filter(!is.na(tpm)) %>% # remove non-detected enzymes for comparison
    select(KO,sampleName,tpm,stoneFormer,pair) %>% 
    compare_counts(count_column = "tpm",
                   grouping_column = "KO") %>% 
    filter(test == "wilcoxon_signed_ranked") %>% 
    mutate(patient = str_split(patient,",") %>% lapply(as.numeric),
           control = str_split(control,",") %>% lapply(as.numeric)) %>% 
    mutate(patient_mean = sapply(patient,function(x){mean(x) %>% round(1)}),
           control_mean = sapply(control,function(x){mean(x) %>% round(1)}),
           patient_se = sapply(patient,function(x){se(x) %>% round(1)}),
           control_se = sapply(control,function(x){se(x) %>% round(1)}),
           patient_label = paste(patient_mean,plusminus,patient_se),
           control_label = paste(control_mean,plusminus,control_se))
}
buta_abundance <- buta_meta %>% compare_abundance() 
ascor_abundance <- ascor_meta %>% compare_abundance()
glyox_abundance <- glyox_meta %>% compare_abundance()
uric_abundance <-uric_meta %>% compare_abundance()
cyst_abundance <- cyst_meta %>% compare_abundance()
  
## 9. add in abundance comparison into presence df
buta_fortable <- buta_presence %>% left_join(buta_abundance %>% select(KO,patient_label,control_label,pval,higher),by="KO")
ascor_fortable <- ascor_presence %>% left_join(ascor_abundance %>% select(KO,patient_label,control_label,pval,higher),by="KO")
glyox_fortable <- glyox_presence %>% left_join(glyox_abundance %>% select(KO,patient_label,control_label,pval,higher),by="KO")
uric_fortable <- uric_presence %>% left_join(uric_abundance %>% select(KO,patient_label,control_label,pval,higher),by="KO")
cyst_fortable <- cyst_presence %>% left_join(cyst_abundance %>% select(KO,patient_label,control_label,pval,higher),by="KO")


## 10. draw table of individual gene's abundance and presence
draw_table <- function(df) {
  df2 <- df %>% 
    arrange(pval) %>% 
    # select(gene_descr,gene_name,control,patient,pval,higher,gene_abbrev,KO,EC) %>% #to check if regex went well
    # select(gene_name,control,patient,pval,higher,gene_abbrev,KO,EC) %>% 
    select(gene_name,n_patient,n_control,patient_label,control_label,pval,higher,KO,gene_abbrev,EC) %>% 
    mutate(pval = round(pval,3)) %>% 
    rename(`Gene name`=gene_name,
           `patient (#)`=n_patient,
           `control (#)`=n_control,
           `patient (tpm)` = patient_label,
           `control (tpm)` = control_label,
           `p value`=pval,
           `Higher in`=higher) %>% 
    as.data.frame()
  
  df2[is.na(df2)] <- "-"
  df2
}
buta_table <- buta_fortable %>% draw_table()
ascor_table <- ascor_fortable %>% draw_table()
glyox_table <- glyox_fortable %>% draw_table()
uric_table <- uric_fortable %>% draw_table()
cyst_table <- cyst_fortable %>% draw_table()

xlsx::write.xlsx(buta_table,"./output/7_three_pathway_genes/buta_table.xlsx",col.names = T,row.names = F)
xlsx::write.xlsx(ascor_table,"./output/7_three_pathway_genes/ascor_table.xlsx",col.names = T,row.names = F)
xlsx::write.xlsx(glyox_table,"./output/7_three_pathway_genes/glyox_table.xlsx",col.names = T,row.names = F)
xlsx::write.xlsx(uric_table,"./output/7_three_pathway_genes/uric_table.xlsx",col.names = T,row.names = F)
xlsx::write.xlsx(cyst_table,"./output/7_three_pathway_genes/cyst_table.xlsx",col.names = T,row.names = F)

### table for thesis
xlsx::write.xlsx(buta_table,"../thesis/chapter1/buta_table_c1.xlsx",col.names = T,row.names = F)
xlsx::write.xlsx(ascor_table,"../thesis/chapter1/ascor_table_c1.xlsx",col.names = T,row.names = F)
xlsx::write.xlsx(glyox_table,"../thesis/chapter1/glyox_table_c1.xlsx",col.names = T,row.names = F)
xlsx::write.xlsx(uric_table,"../thesis/chapter1/uric_table_c1.xlsx",col.names = T,row.names = F)
xlsx::write.xlsx(cyst_table,"../thesis/chapter1/cyst_table_c1.xlsx",col.names = T,row.names = F)

# Full pathway comparison

## 11. Calculate percentage of abundance for an entire pathway (out of a total of 1000000 tpm per sample)
pwy_total <- function(df) {
  df %>% 
    select(pathway_name,KO,sampleName,stoneFormer,pair,tpm) %>% 
    group_by(sampleName,stoneFormer,pair) %>% 
    summarize(sum_tpm = sum(tpm,na.rm = T),
              pathway_name = unique(pathway_name)) %>% 
    mutate(perc_tpm = sum_tpm / 1000000 * 100) %>% 
    mutate(pathway_label = pathway_name %>% str_replace("\\[.*\\]",""))
}

buta_pwy_total <- buta_meta %>% pwy_total()
glyox_pwy_total <- glyox_meta %>% pwy_total()  
ascor_pwy_total <- ascor_meta %>% pwy_total()  
uric_pwy_total <- uric_meta %>% pwy_total()
cyst_pwy_total <- cyst_meta %>% pwy_total()

(five_pwy_total_hist <- rbind(buta_pwy_total,
                               ascor_pwy_total,
                               glyox_pwy_total,
                              uric_pwy_total,
                              cyst_pwy_total) %>% 
    ggplot(aes(x=perc_tpm)) + geom_histogram(bins=8) + 
    facet_wrap(~pathway_name,scales="free_x")) # distribution looks mostly normal


## 12. Compare percentage of total pathway abundance (% total tpm)
buta_pwy_total_compared <- buta_pwy_total %>% compare_counts(count_column = "perc_tpm",grouping_column = "pathway_label") 
glyox_pwy_total_compared <- glyox_pwy_total %>% compare_counts(count_column = "perc_tpm",grouping_column = "pathway_label") 
ascor_pwy_total_compared <- ascor_pwy_total %>% compare_counts(count_column = "perc_tpm",grouping_column = "pathway_label") 
uric_pwy_total_compared <- uric_pwy_total %>% compare_counts(count_column = "perc_tpm",grouping_column = "pathway_label") 
cyst_pwy_total_compared <- cyst_pwy_total %>% compare_counts(count_column = "perc_tpm",grouping_column = "pathway_label") 


five_pwy_total_compared <- rbind(buta_pwy_total_compared,
                                  glyox_pwy_total_compared,
                                  ascor_pwy_total_compared,
                                 uric_pwy_total_compared,
                                 cyst_pwy_total_compared) %>% 
  filter(test == "ttest") # we select paired t test because the distribution looks normal in 11.


## table for thesis
five_pwy_total_compared_table <- five_pwy_total_compared %>% 
  ungroup() %>% 
  mutate(pval = round(pval,3)) %>% 
  mutate(patient = str_split(patient,",") %>% lapply(as.numeric),
         control = str_split(control,",") %>% lapply(as.numeric)) %>% 
  mutate(patient_mean = sapply(patient,function(x){mean(x) %>% round(3)}),
         control_mean = sapply(control,function(x){mean(x) %>% round(3)}),
         patient_se = sapply(patient,function(x){se(x) %>% round(3)}),
         control_se = sapply(control,function(x){se(x) %>% round(3)})) %>% 
  mutate(patient_label = paste(patient_mean,"\u00b1",patient_se),
         control_label = paste(control_mean,"\u00b1",control_se)) %>% 
  select(pathway_label,patient_label,control_label,pval,higher) %>% 
  rename(`Metabolic Pathway` = pathway_label,
         patient = patient_label,
         control = control_label,
         `p value` = pval,
         `higher in` = higher) %>% 
  as.data.frame()

xlsx::write.xlsx(five_pwy_total_compared_table,file = "../thesis/chapter1/five_pwy_total_compared_table_c1.xlsx",col.names = T,row.names = F)
xlsx::write.xlsx(five_pwy_total_compared_table,file = "./output/7_three_pathway_genes/five_pwy_total_compared_table.xlsx",col.names = T,row.names = F)


## 13. Plot percentage of total pathway abundance (% total tpm)
plot_paired <- function(df,value,boxplot=F,facet_scales) {
  base_plot <- df %>% 
    group_by(pair) %>% 
    
    mutate_("values" = value) %>% 
    mutate(compare_values = paste0(values,collapse=","),
           compare_grp = paste0(stoneFormer,collapse=","),
           compare_pair = paste0(pair,collapse=",")) %>% 
    mutate(higher = mapply(which_higher,
                           comma_separated_num = compare_values,
                           comma_separated_grp = compare_grp)) %>% 
    ungroup() %>% 
    ggplot(aes(x=stoneFormer,y=values)) 
  
  if(boxplot) {
    base_plot +
      geom_boxplot()+
      facet_grid(~pathway_name)
  } else {
    base_plot +
      geom_point(aes(color=stoneFormer)) +
      geom_line(aes(group = pair,color=higher,alpha=0.2)) +
      scale_color_manual(values=c("#00BFC4","#F8766D","black")) +
      scale_alpha_continuous(guide=F)+
      facet_wrap("pathway_name",scales=facet_scales)
  }
}

(pwy_total_plot <- rbind(buta_pwy_total,
                         ascor_pwy_total,
                         glyox_pwy_total,
                         uric_pwy_total,
                         cyst_pwy_total) %>% 
    mutate(pathway_name = str_wrap(pathway_name,20)) %>% 
    plot_paired(boxplot = F,value="perc_tpm",facet_scales = "fixed") +
    
    theme1 +
    ylab("Percentage of total TPM ( % TPM )") +
    xlab("") +
    labs(color="") +
    theme(plot.margin = margin(0,10,0,0),
          strip.text = element_text(size=6),
          axis.title.y = element_text(size=8),
          legend.position = "none")) 

ggsave(filename = "./output/7_three_pathway_genes/pwy_total_plot.tiff",
       pwy_total_plot,
       dpi = 280,
       width = plos_dims(10,8,0.7)$plos.width,
       height = plos_dims(10,8,0.7)$plos.height,
       units = "cm")

### 14. Calculate percentage of presence for entire pathway
pwy_presence <- function(df) {
  df %>% mutate(present = (tpm > 0 & !is.na(tpm)) %>% as.numeric()) %>% 
    group_by(pathway_name,pair,stoneFormer,sampleName) %>% 
    mutate(total_pwy_genes = n()) %>% 
    summarize(sum_present = sum(present),
              total_pwy_genes = unique(total_pwy_genes),
              percent_coverage = (sum_present/total_pwy_genes * 100) %>% round(1)) %>% 
    mutate(pathway_label = pathway_name %>% str_replace("\\[.*\\]","") %>% paste("(",total_pwy_genes,"genes total )")) %>% 
    ungroup()
}

buta_pwy_presence <- buta_meta %>% pwy_presence()
glyox_pwy_presence <- glyox_meta %>% pwy_presence() 
ascor_pwy_presence <- ascor_meta %>% pwy_presence()
uric_pwy_presence <- uric_meta %>% pwy_presence()
cyst_pwy_presence <- cyst_meta %>% pwy_presence()

### 15. compare percentage of presence 
buta_pwy_presence_compared <- buta_pwy_presence %>% compare_counts(count_column = "percent_coverage",grouping_column = "pathway_label")
glyox_pwy_presence_compared <- glyox_pwy_presence %>% compare_counts(count_column = "percent_coverage",grouping_column = "pathway_label")
ascor_pwy_presence_compared <- ascor_pwy_presence %>% compare_counts(count_column = "percent_coverage",grouping_column = "pathway_label")
uric_pwy_presence_compared <- uric_pwy_presence %>% compare_counts(count_column = "percent_coverage",grouping_column = "pathway_label")
cyst_pwy_presence_compared <- cyst_pwy_presence %>% compare_counts(count_column = "percent_coverage",grouping_column = "pathway_label")

five_pwy_presence_compared <- rbind(buta_pwy_presence_compared,
                                  glyox_pwy_presence_compared,
                                  ascor_pwy_presence_compared,
                                  uric_pwy_presence_compared,
                                  cyst_pwy_presence_compared) %>% 
  filter(test == "ttest") # we select paired t test because the presence value distribution looks normal.

## table for thesis
five_pwy_presence_compared_table <- five_pwy_presence_compared %>% 
  ungroup() %>% 
  mutate(pval = round(pval,3)) %>% 
  mutate(patient = str_split(patient,",") %>% lapply(as.numeric),
         control = str_split(control,",") %>% lapply(as.numeric)) %>% 
  mutate(patient_mean = sapply(patient,function(x){mean(x) %>% round(1)}),
         control_mean = sapply(control,function(x){mean(x) %>% round(1)}),
         patient_se = sapply(patient,function(x){se(x) %>% round(1)}),
         control_se = sapply(control,function(x){se(x) %>% round(1)})) %>% 
  mutate(patient_label = paste(patient_mean,"\u00b1",patient_se),
         control_label = paste(control_mean,"\u00b1",control_se)) %>% 
  select(pathway_label,patient_label,control_label,pval,higher) %>% 
  rename(`Metabolic Pathway` = pathway_label,
         patient = patient_label,
         control = control_label,
         `p value` = pval,
         `higher in` = higher) %>% 
  as.data.frame()

xlsx::write.xlsx(five_pwy_presence_compared_table,file = "../thesis/chapter1/five_pwy_presence_compared_table_c1.xlsx",col.names = T,row.names = F)

xlsx::write.xlsx(five_pwy_presence_compared_table,file = "./output/7_three_pathway_genes/five_pwy_presence_compared_table.xlsx",col.names = T,row.names = F)

### 16. plot percentage of presence
(pwy_presence_plot <- rbind(buta_pwy_presence,
                            glyox_pwy_presence,
                            ascor_pwy_presence,
                            uric_pwy_presence,
                            cyst_pwy_presence) %>% 
    mutate(pathway_name = str_wrap(pathway_name,20)) %>%
    mutate(pathway_name = paste0(pathway_name,"\n",total_pwy_genes," genes in pathway")) %>% 
    plot_paired(value = "sum_present",boxplot = F,facet_scales = "free_y") +
    theme1 +
    ylab("number of genes detected") +
    xlab("") +
    labs(color="") +
    theme(plot.margin = margin(0,10,0,0),
          strip.text = element_text(size=5),
          axis.title = element_text(size=6),
          axis.text = element_text(size=6),
          legend.position = "none"))

ggsave(filename = "./output/7_three_pathway_genes/pwy_presence_plot.tiff",
       pwy_presence_plot,
       dpi = 280,
       width = 18,
       height = 9,
       units = "cm")





# check all
all_tpm <- kegg_tpm_results
all_meta <- all_tpm %>% add_metadata()
# all_abundance <- all_meta %>% compare_abundance()
# write_rds(all_abundance,"./output/7_three_pathway_genes/all_abundance.rds")
all_abundance <- read_rds("./output/7_three_pathway_genes/all_abundance.rds")

all_hierarchized <- all_abundance %>% 
  left_join(kegg_hierarchy_ec,by="KO") %>% 
  select(-hierarchy1,-hierarchy2,-pathway) %>% 
  unique() %>% 
  group_by_at(vars(-pathway_name)) %>% 
  summarize(pathway_name = paste(pathway_name,collapse=","))

# note: nope, no interesting trend


## Sanity check: How should percentage of total tpm be defined?
## Should we take into account predicted ORFs with no annotations?
## temporary choice: define total tpm as all tpm (with and without annotations)
## what is tpm? http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

orf_stats_results <- read.delim("~/aNALYSIS_PIPELINE/data/metagenomics/MC3_prep_rpkmoutput/orf_stats_results.txt",header = T,sep = "\t")

orf_stats_meta <- orf_stats_results %>% 
  mutate(sampleName = toupper(sampleName)) %>% 
  inner_join(metadata_basic,by = "sampleName")

(orf_stats_overview <- orf_stats_meta %>% 
  ggplot(aes(x=sampleName,y=perc_tpm,fill=status)) +
  geom_bar(stat="identity"))

orf_stats_compared <- orf_stats_meta %>% 
  compare_counts(count_column = "perc_tpm",
                 grouping_column = "status") 

(orf_stats_plot <-
  orf_stats_meta %>% 
  group_by(status,pair) %>% 
  mutate(compare_tpm = paste0(perc_tpm,collapse=","),
         compare_grp = paste0(stoneFormer,collapse=","),
         compare_pair = paste0(pair,collapse=",")) %>% 
  mutate(higher = mapply(which_higher,
                         comma_separated_num = compare_tpm,
                         comma_separated_grp = compare_grp)) %>%
  ungroup() %>%
  ggplot(aes(x=stoneFormer,y=perc_tpm)) +
  geom_point(aes(color=stoneFormer)) +
  geom_line(aes(group = pair,color=higher)) +
  scale_color_manual(values=c("#00BFC4","#F8766D","black")) +
  facet_wrap(~status,nrow=1))



# ###OLD
# #color scheme
# library(scales)
# show_col(hue_pal()(3))
# show_col(hue_pal()(2))
# 
# (kegg_plot_buta <- kegg_rdyplot_buta %>% 
#     filter(dom_group=="control") %>%
#     filter(!str_detect(gene_descr,"por[ABDG]{1}")) %>% 
#     filter(!str_detect(gene_descr,"succinate")) %>% 
#     mutate(gene_label = str_wrap(gene_label,width = 25)) %>% 
#     mutate(stoneFormer = factor(stoneFormer,levels=c("control","patient"))) %>%
#     # ggplot(aes(x=stoneFormer,y=log(ko_tpm + 1),color=stoneFormer)) + 
#     ggplot(aes(x=stoneFormer,y=ko_tpm,color=stoneFormer)) + 
#     geom_boxplot() +
#     geom_point() +
#     facet_wrap(~gene_label,scales="free_y",nrow = 1)+
#     # ggtitle("PT/CTRL differences in the Butanoate pathway (KEGG)") +
#     ggtitle("") +
#     xlab("") +
#     ylab("Gene abundance ( TPM )")+
#     theme_light() +
#     theme_tufte() +
#     theme(legend.position="none") +
#     theme(strip.text = element_text(size=12),
#           axis.title.x = element_text(size = 15,margin = margin(0,0,0,0)),
#           axis.title.y = element_text(size = 15,margin = margin(0,30,0,0)),
#           # axis.text.x = element_text(size = 10,angle=-45,vjust=1,hjust=0),
#           # axis.text.x = element_text(size = 10),
#           axis.text.x = element_text(size = 10),
#           axis.text.y = element_text(size = 15),
#           axis.line = element_line(colour = "black"),
#           panel.spacing = unit(3,"lines"),
#           plot.margin = unit(c(3,5,5,5),"lines")) +
#     # theme(strip.background =element_rect(fill="white"))+
#     theme(strip.text.x = element_text(colour = 'black',margin = margin(0,0,20,0))) +
#     scale_color_manual(values=c("#00BFC4","#F8766D"))
# )
# 
# pdf("~/aNALYSIS_PIPELINE/git_world/gut_metagenomics/publication/3_butanoate_enzymes(done)/butanoate_enzymes.pdf",width = 10,height=7)
# kegg_plot_buta
# dev.off()
# 
# 
# (kegg_plot_pfor <- kegg_rdyplot_buta %>% 
#   filter(dom_group=="control") %>%
#   filter(str_detect(gene_descr,"por[ABDG]")) %>% 
#     mutate(gene_label = str_wrap(gene_label,width = 19)) %>% 
#     mutate(stoneFormer = factor(stoneFormer,levels=c("control","patient"))) %>%
#     # ggplot(aes(x=stoneFormer,y=log(ko_tpm + 1),color=stoneFormer)) + 
#     ggplot(aes(x=stoneFormer,y=ko_tpm,color=stoneFormer)) + 
#     geom_boxplot() +
#     geom_point() +
#     facet_wrap(~gene_label,scales="free_y",nrow = 1)+
#     # ggtitle("PT/CTRL differences in the Butanoate pathway (KEGG)") +
#     ggtitle("") +
#     xlab("") +
#     ylab("Gene abundance ( TPM )")+
#     theme_light() +
#     theme_tufte() +
#     theme(legend.position="none") +
#     theme(strip.text = element_text(size=12),
#           axis.title.x = element_text(size = 15,margin = margin(0,0,0,0)),
#           axis.title.y = element_text(size = 15,margin = margin(0,30,0,0)),
#           # axis.text.x = element_text(size = 10,angle=-45,vjust=1,hjust=0),
#           # axis.text.x = element_text(size = 10),
#           axis.text.x = element_text(size = 10),
#           axis.text.y = element_text(size = 15),
#           axis.line = element_line(colour = "black"),
#           panel.spacing = unit(3,"lines"),
#           plot.margin = unit(c(3,5,5,5),"lines")) +
#     # theme(strip.background =element_rect(fill="white"))+
#     theme(strip.text.x = element_text(colour = 'black',margin = margin(0,0,20,0)))+
#     scale_color_manual(values=c("#00BFC4","#F8766D"))
# )
# 
# pdf("~/aNALYSIS_PIPELINE/git_world/gut_metagenomics/publication/3_butanoate_enzymes(done)/pfor_subunits.pdf",width = 10,height=7)
# kegg_plot_pfor
# dev.off()
# 
# #color scheme
# library(scales)
# show_col(hue_pal()(3))
# show_col(hue_pal()(2))
# #PAIRED PLOTS
# 
# # butanoate enzymes higher in ctrl
# (kegg_pairedplot_buta_ctrl <- kegg_wilcox_buta %>% filter(pval < 0.1) %>%
#     mutate(pair=rep(list(1:17),nrow(.))) %>% 
#     unnest(pair,patient,control) %>% 
#     mutate(higher = ifelse(patient > control,"patient",
#                            ifelse(patient == control,"neutral","control"))) %>%
#     mutate(higher = factor(higher,levels=c("control","patient","neutral"))) %>%
#     gather("stoneFormer","TPM",7:8) %>% 
#     filter(dom_group == "control") %>% 
#     filter(!str_detect(gene_descr,"por")) %>% 
#     filter(!str_detect(gene_descr,"succinate")) %>% 
#     mutate(gene_label = str_replace(gene_descr,".*[aA-zZ]+;","")) %>% 
#     mutate(gene_label = str_replace(gene_label,"\\[.*\\]","")) %>% 
#     mutate(gene_label = str_wrap(gene_label,width=30)) %>% 
#     ggplot(aes(x=stoneFormer,y=TPM,group=pair,color=higher)) +
#     geom_point() +
#     geom_line() +
#     facet_wrap(~gene_label, scales = "free_y") +
#     ggtitle("") +
#     xlab("") +
#     ylab("Gene abundance ( TPM )")+
#     theme_light() +
#     theme_tufte() +
#     theme(legend.position="none") +
#     theme(strip.text = element_text(size=12),
#           axis.title.x = element_text(size = 15,margin = margin(0,0,0,0)),
#           axis.title.y = element_text(size = 15,margin = margin(0,30,0,0)),
#           # axis.text.x = element_text(size = 10,angle=-45,vjust=1,hjust=0),
#           # axis.text.x = element_text(size = 10),
#           axis.text.x = element_text(size = 10),
#           axis.text.y = element_text(size = 15),
#           axis.line = element_line(colour = "black"),
#           panel.spacing = unit(3,"lines"),
#           plot.margin = unit(c(3,5,5,5),"lines")) +
#     # theme(strip.background =element_rect(fill="white"))+
#     theme(strip.text.x = element_text(colour = 'black',margin = margin(0,0,20,0))) +
#     scale_color_manual(values=c("#00BFC4","#F8766D", "black"))
# )
# 
# pdf("~/aNALYSIS_PIPELINE/git_world/gut_metagenomics/publication/3_butanoate_enzymes(done)/paired_butanoate_enzymes_ctrl.pdf",width = 12,height=7)
# kegg_pairedplot_buta_ctrl
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# # PFOR enzymes higher in ctrl
# (kegg_pairedplot_pfor <- kegg_wilcox_buta %>% 
#     # filter(pval < 0.1) %>%
#     filter(str_detect(gene_descr,"por")) %>% 
#     mutate(pair=rep(list(1:17),nrow(.))) %>% 
#     unnest(pair,patient,control) %>% 
#     mutate(higher = ifelse(patient > control,"patient",
#                            ifelse(patient == control,"neutral","control"))) %>%
#     mutate(higher = factor(higher,levels=c("control","patient","neutral"))) %>%
#     gather("stoneFormer","TPM",7:8) %>% 
#     filter(dom_group == "control") %>% 
#     filter(str_detect(gene_descr,"por")) %>% 
#     mutate(gene_label = str_replace(gene_descr,".*[aA-zZ]+;","")) %>% 
#     mutate(gene_label = str_replace(gene_label,"\\[.*\\]","")) %>% 
#     mutate(gene_label = str_wrap(gene_label,width=19)) %>% 
#     ggplot(aes(x=stoneFormer,y=TPM,group=pair,color=higher)) +
#     geom_point() +
#     geom_line() +
#     facet_wrap(~gene_label, scales = "free_y", nrow=1) +
#     ggtitle("") +
#     xlab("") +
#     ylab("Gene abundance ( TPM )")+
#     theme_light() +
#     theme_tufte() +
#     theme(legend.position="none") +
#     theme(strip.text = element_text(size=12),
#           axis.title.x = element_text(size = 15,margin = margin(0,0,0,0)),
#           axis.title.y = element_text(size = 15,margin = margin(0,30,0,0)),
#           # axis.text.x = element_text(size = 10,angle=-45,vjust=1,hjust=0),
#           # axis.text.x = element_text(size = 10),
#           axis.text.x = element_text(size = 10),
#           axis.text.y = element_text(size = 15),
#           axis.line = element_line(colour = "black"),
#           panel.spacing = unit(3,"lines"),
#           plot.margin = unit(c(3,3,3,3),"lines")) +
#     # theme(strip.background =element_rect(fill="white"))+
#     theme(strip.text.x = element_text(colour = 'black',margin = margin(0,0,20,0))) +
#     scale_color_manual(values=c("#00BFC4","#F8766D", "black"))
# )
# 
# pdf("~/aNALYSIS_PIPELINE/git_world/gut_metagenomics/publication/3_butanoate_enzymes(done)/paired_pfor_enzymes.pdf",width = 12,height=7)
# kegg_pairedplot_pfor
# dev.off()
# 
# # butanoate enzymes higher in pt
# (kegg_pairedplot_buta_pt <- kegg_wilcox_buta %>% filter(pval < 0.1) %>%
#     mutate(pair=rep(list(1:17),nrow(.))) %>% 
#     unnest(pair,patient,control) %>% 
#     mutate(higher = ifelse(patient > control,"patient",
#                            ifelse(patient == control,"neutral","control"))) %>%
#     mutate(higher = factor(higher,levels=c("control","patient","neutral"))) %>%
#     gather("stoneFormer","TPM",7:8) %>% 
#     filter(dom_group == "patient") %>% 
#     mutate(gene_label = str_replace(gene_descr,".*[aA-zZ]+;","")) %>% 
#     mutate(gene_label = str_replace(gene_label,"\\[.*\\]","")) %>% 
#     mutate(gene_label = str_wrap(gene_label,width=30)) %>% 
#     ggplot(aes(x=stoneFormer,y=TPM,group=pair,color=higher)) +
#     geom_point()+
#     geom_line() +
#     facet_wrap(~gene_label, scales = "free_y") +
#     ggtitle("") +
#     xlab("") +
#     ylab("Gene abundance ( TPM )")+
#     theme_light() +
#     theme_tufte() +
#     theme(legend.position="none") +
#     theme(strip.text = element_text(size=12),
#           axis.title.x = element_text(size = 15,margin = margin(0,0,0,0)),
#           axis.title.y = element_text(size = 15,margin = margin(0,30,0,0)),
#           # axis.text.x = element_text(size = 10,angle=-45,vjust=1,hjust=0),
#           # axis.text.x = element_text(size = 10),
#           axis.text.x = element_text(size = 10),
#           axis.text.y = element_text(size = 15),
#           axis.line = element_line(colour = "black"),
#           panel.spacing = unit(3,"lines"),
#           plot.margin = unit(c(3,5,5,5),"lines")) +
#     # theme(strip.background =element_rect(fill="white"))+
#     theme(strip.text.x = element_text(colour = 'black',margin = margin(0,0,20,0))) +
#     scale_color_manual(values=c("#00BFC4","#F8766D", "black"))
# )
# 
# pdf("~/aNALYSIS_PIPELINE/git_world/gut_metagenomics/publication/3_butanoate_enzymes(done)/paired_butanoate_enzymes_pt.pdf",width = 12,height=7)
# kegg_pairedplot_buta_pt
# dev.off()
# 
# 
# (kegg_pairedplot_glyox <- kegg_wilcox_glyox %>% filter(pval < 0.1) %>%
#     mutate(pair=rep(list(1:17),nrow(.))) %>% 
#     unnest(pair,patient,control) %>% 
#     mutate(higher = ifelse(patient > control,"patient",
#                            ifelse(patient == control,"neutral","control"))) %>%
#     mutate(higher = factor(higher,levels=c("control","patient","neutral"))) %>%
#     gather("stoneFormer","TPM",7:8) %>% 
#     # filter(dom_group == "patient") %>% 
#     mutate(gene_label = str_replace(gene_descr,".*[aA-zZ]+;","")) %>% 
#     mutate(gene_label = str_replace(gene_label,"E[0-9]{1,3}\\.[0-9]{1,3}\\.[0-9]{1,3}\\.[0-9]{1,3};","")) %>% 
#     mutate(gene_label = str_replace(gene_label,"\\[.*\\]","")) %>% 
#     mutate(gene_label = str_wrap(gene_label,width=30)) %>% 
#     ggplot(aes(x=stoneFormer,y=TPM,group=pair,color=higher)) +
#     geom_point()+
#     geom_line() +
#     facet_wrap(~gene_label, scales = "free_y",nrow=1) +
#     ggtitle("") +
#     xlab("") +
#     ylab("Gene abundance ( TPM )")+
#     theme_light() +
#     theme_tufte() +
#     theme(legend.position="none") +
#     theme(strip.text = element_text(size=12),
#           axis.title.x = element_text(size = 15,margin = margin(0,0,0,0)),
#           axis.title.y = element_text(size = 15,margin = margin(0,30,0,0)),
#           # axis.text.x = element_text(size = 10,angle=-45,vjust=1,hjust=0),
#           # axis.text.x = element_text(size = 10),
#           axis.text.x = element_text(size = 10),
#           axis.text.y = element_text(size = 15),
#           axis.line = element_line(colour = "black"),
#           panel.spacing = unit(3,"lines"),
#           plot.margin = unit(c(3,5,5,5),"lines")) +
#     # theme(strip.background =element_rect(fill="white"))+
#     theme(strip.text.x = element_text(colour = 'black',margin = margin(0,0,20,0))) +
#     scale_color_manual(values=c("#00BFC4","#F8766D", "black"))
# )
# 
# pdf("~/aNALYSIS_PIPELINE/git_world/gut_metagenomics/publication/4_ascorbate_glyxoylate_enzymes(done)//paired_glyox_enzymes.pdf",width = 15,height=7)
# kegg_pairedplot_glyox
# dev.off()
# 
# 
# 
# # # Compare overall pathway
# # write_rds(kegg_tpm_overview,"~/aNALYSIS_PIPELINE/data/metagenomics/MC3_prep_rpkmoutput/kegg_tpm_overview.rds")
# # write_rds(kegg_tpm_results_meta,"~/aNALYSIS_PIPELINE/data/metagenomics/MC3_prep_rpkmoutput/kegg_tpm_results_meta.rds")
# 
# 
# # EXTRAS: Draw table for high/low to check consistency of results across all pt/ctrl pairs
# factor_order1 <- c(c(4,16),c(2,1,5,14),c(7,9,11,17),c(6,10,12,13),c(15),c(3,8))
# 
# kegg_tpm_meta %>% 
#   ungroup() %>% 
#   left_join(kegg_hierarchy,by="KO") %>% 
#   filter(str_detect(pathway_name,"Butanoate")) %>% 
#   filter(KO %in% ko_candidates_buta) %>% 
#   select(stoneFormer,pair,gene_descr,ko_tpm) %>%
#   spread(stoneFormer,ko_tpm) %>% 
#   mutate(higher = ifelse(patient > control," ","control")) %>% 
#   select(-patient,-control) %>% 
#   mutate(pair = factor(pair,levels = factor_order1)) %>% 
#   # mutate(pair = as.integer(pair)) %>% 
#   group_by(gene_descr) %>% 
#   mutate(sum_ctrl = sum(higher=="control")) %>% 
#   ungroup() %>% 
#   filter(sum_ctrl > 10) %>% 
#   mutate(gene_descr = reorder(gene_descr,sum_ctrl)) %>% 
#   spread(pair,higher) %>% 
#   View()
# 
