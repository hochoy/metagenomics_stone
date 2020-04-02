# Combine protein sequences for the 3 subunits of Oxalate Oxidoreductase into a single fasta file and build it into a database file.

## 0. set working directory, set path, make output directory
setwd("~/Analysis_of_gut/metagenomics/")
Sys.setenv(PATH = readLines("scripts_metagenomics/0_setup/bash_path.txt"))
system("mkdir ~/Analysis_of_gut/metagenomics/output/6_blast_for_oxalate_oxidoreductase/")
system("mkdir ~/Analysis_of_gut/metagenomics/output/6_blast_for_oxalate_oxidoreductase/diamond_db")
system("mkdir ~/Analysis_of_gut/metagenomics/output/6_blast_for_oxalate_oxidoreductase/diamond_alignment")

## 1. load libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(R.utils)
source("./scripts_metagenomics/0_setup/0_functions.R")

## 2.Concatenate the protein.fasta files
cat_fasta <- function(fasta_dir,output_file) {
    cmd <- paste0("cat ",fasta_dir,"*.fasta"," > ",output_file)
    system(cmd)
}
# cat_fasta(fasta_dir = "./input/blast/oxalate_oxidoreductase/",
#           output_file = "./output/6_blast_for_oxalate_oxidoreductase/diamond_db/oxalateoxidoreductase_subunits.fasta")


## 3. linearize protein fasta file (protein sequences are currently wrapped: meaning there's a newline separating each protein seq line)
linearize_fasta <- function(fasta_file,output_file) {
  bash_script <- "./scripts_metagenomics/6_blast_for_oxalate_oxidoreductase/bash_scripts/2_linearize_pfasta.sh"
  cmd <- paste(bash_script,
               "-i",fasta_file,
               "-o",output_file,
               "-r")
  system(cmd)
}
# linearize_fasta(fasta_file = "./output/6_blast_for_oxalate_oxidoreductase/diamond_db/oxalateoxidoreductase_subunits.fasta",
#                 output_file = "./output/6_blast_for_oxalate_oxidoreductase/diamond_db/oxalateoxidoreductase_linear.fasta")


## 4. deduplicate protein sequences (two or more sequences that have identical headers). not really required in this case but just leave it in.
deduplicate_fasta <- function(linear_fasta,output_file) {
  bash_script <- "./scripts_metagenomics/6_blast_for_oxalate_oxidoreductase/bash_scripts/3_unique_pfasta.sh"
  cmd <- paste(bash_script,
               "-i", linear_fasta,
               "-o", output_file,
               "-d", dirname(linear_fasta),
               "-r")
  system(cmd)
  
}
# deduplicate_fasta(linear_fasta = "./output/6_blast_for_oxalate_oxidoreductase/diamond_db/oxalateoxidoreductase_linear.fasta",
#                   output_file = "./output/6_blast_for_oxalate_oxidoreductase/diamond_db/oxalateoxidoreductase_unique.fasta")


## 5. grep protein fasta headers
headers_fasta <- function(fasta_file,output_file){
  bash_script <- "./scripts_metagenomics/6_blast_for_oxalate_oxidoreductase/bash_scripts/4_header_pfasta.sh"
  cmd <- paste(bash_script,
               "-i", fasta_file,
               "-o", output_file,
               "-r")
  system(cmd)
}
# headers_fasta(fasta_file = "./output/6_blast_for_oxalate_oxidoreductase/diamond_db/oxalateoxidoreductase_unique.fasta",
#               output_file = "./output/6_blast_for_oxalate_oxidoreductase/diamond_db/oxalateoxidoreductase_header.fasta")



# 5. makedb using protein fasta files
makedb_fasta <-function(fasta_file,output_file) {
  bash_script <- "./scripts_metagenomics/6_blast_for_oxalate_oxidoreductase/bash_scripts/5_makedb_pfasta.sh"
  cmd <- paste(bash_script,
               "-i",fasta_file,
               "-o",output_file,
               "-r")
  system(cmd)
  
}
# makedb_fasta(fasta_file = "./output/6_blast_for_oxalate_oxidoreductase/diamond_db/oxalateoxidoreductase_unique.fasta",
#              output_file = "./output/6_blast_for_oxalate_oxidoreductase/diamond_db/oxalateoxidoreductase_unique")


# 6. run diamond's blastx
run_diamond <- function(db, dna_file, output_file, setting="default") {
  bash_script <- "./scripts_metagenomics/6_blast_for_oxalate_oxidoreductase/bash_scripts/6_blastx_dfasta.sh"
  cmd <- paste(bash_script,
               "-i",dna_file,
               "-d", db,
               "-o", output_file,
               "-s", setting,
               "-r")
  system(cmd)
}
#test run on single sample
# run_diamond(dna_file = "./input/clean_reads/CTRLXXX_merged.gz",
#             db = "./output/6_blast_for_oxalate_oxidoreductase/diamond_db/oxalateoxidoreductase_unique",
#             output_file = "./output/6_blast_for_oxalate_oxidoreductase/diamond_alignment/alignment_CTRLXXX_default",
#             setting = "default")


# 7. build dataframe to organize all samples
reads_folder <- "./input/clean_reads/"
reads_files <- list.files(reads_folder,"*.gz",full.names = T)

diamond_df <- data.frame(reads_files = reads_files) %>% 
  mutate(sampleName = str_extract(reads_files,"CTRLA[0-9,A-Z]+|PT[0-9,A-Z]+") %>% trimws(),
         alignment_files = paste0("./output/6_blast_for_oxalate_oxidoreductase/diamond_alignment/",
                                  sampleName,
                                  "_alignment_default"))

# 8. Check for already completed alignments. 
done_filenames <-list.files(path = "./output/6_blast_for_oxalate_oxidoreductase/diamond_alignment/",full.names = T)
done_sizes <- file.info(done_filenames)$size
done_samples <- data.frame(done_filenames = done_filenames,
                               done_sizes = done_sizes) %>% 
  mutate(alignment_complete = done_sizes > 1000) %>% 
  filter(alignment_complete) %>% 
  mutate(sampleName = str_extract(done_filenames,"CTRLA[0-9,A-Z]{1,3}|PT[0-9,A-Z]{1,3}") %>% trimws())

# 9. Filter out already completed alignments.
diamond_df2 <- diamond_df %>% 
  filter(!sampleName %in% done_samples$sampleName)
nrow(diamond_df) == (nrow(done_samples) + nrow(diamond_df2))

# 10. run all unaligned samples
diamond_df2 %>% 
  mutate(run = mapply(FUN=run_diamond,
                      dna_file = reads_files,
                      db = "./output/6_blast_for_oxalate_oxidoreductase/diamond_db/oxalateoxidoreductase_unique",
                      output_file = alignment_files))




