# Prepare metagenomics metadata

## 0. Set working directory, make output directory
setwd("~/Analysis_of_gut/metagenomics/")
bash_path <- readLines("scripts_metagenomics/0_setup/bash_path.txt")
Sys.setenv(PATH = bash_path)
system("mkdir output/1_metadata")
system("mkdir ../thesis/")

## 1. Load required libraries
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(xlsx)
source("./scripts_metagenomics/0_setup/0_functions.R")


## 2. First load in the participant metadata and additional 16s and metagenomics metadata
metadata <- read.csv("~/Analysis_of_gut/metagenomics/input/metadata/ptctrl_metadata.csv")
metadata_momics <- read.csv("~/Analysis_of_gut/metagenomics/input/metadata/metagenomics_metadata.csv")
metadata_stonetype <- read.xlsx("~/Analysis_of_gut/metagenomics_fixed/input/metadata/Stone Composition Information.xlsx",1,header = T)
# stonetype_pair <- metadata %>% 
#   select(sampleName,pair) %>% 
#   left_join(metadata_stonetype,by="sampleName") %>%
#   select(pair,Stone.Compositon) %>% 
#   filter(str_detect(Stone.Compositon,"CaOx|uric"))

## 3. A single woodrat fecal sample was included in the metadata. Do not include it.
##  Filter out the woodrat entry from metadata
metadata <- metadata %>% filter(sampleName != "WOODRAT") %>% droplevels()
metadata_momics <- metadata_momics %>% filter(sampleName != "WOODRAT") %>% droplevels()
metadata_final <- metadata %>% 
  left_join(metadata_momics,by="sampleName") #%>% 
  # filter(pair %in% stonetype_pair$pair)

write.csv(metadata_final,
          "./output/1_metadata/metadata_final.csv",
          row.names = F)

## 4. summarize metadata
se <- function(x) {sqrt(var(x)/length(x)) %>% round(1)}

gender_df <- metadata_final %>% 
  select(stoneFormer,gender) %>% 
  group_by(stoneFormer,gender) %>% 
  summarize(count=n()) %>% 
  rename(metric = gender) %>% 
  spread(stoneFormer,count) %>% 
  arrange(desc(metric))

age_df <- metadata_final %>% 
  select(stoneFormer,age) %>% 
  group_by(stoneFormer) %>% 
  summarize(age_mean = mean(age) %>% round(0),
            age_se = se(age),
            age_label = paste(age_mean,"\u00b1",age_se)) %>% 
  mutate(metric = "Mean age") %>% 
  select(stoneFormer,age_label,metric) %>% 
  spread(stoneFormer,age_label)
  
metadata_table <- rbind(gender_df,
                        age_df) %>% 
  rename(` `=metric,
         Patient = patient,
         Control = control) %>% 
  select(` `,Patient,Control) %>% 
  as.data.frame() 

write.xlsx(metadata_table,"../thesis/chapter1/metadata_table_c1.xlsx",col.names = T,row.names = F)
