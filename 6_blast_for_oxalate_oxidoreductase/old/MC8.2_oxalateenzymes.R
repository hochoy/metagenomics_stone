library(R.utils)
#tidyverse
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(reshape2)
library(plotly)



# read in mapping file
oxalate_db_header <- read_lines("~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_header.fasta") 

oxalate_db_map <- data.frame(header = oxalate_db_header) %>% 
  mutate(header = as.character(header),
         prot_id = str_extract(header,">[^ ]+ ") %>% trimws(),
         description = str_replace(header,">[^ ]+ ","") %>% trimws(),
         nchar_id = nchar(prot_id)) %>% 
  mutate(subject = prot_id %>% str_replace(">",""))



# read in TEST diamond blastx output

CTRLXXX_oxalateproteins_default_raw <- read_delim("~/custom_diamond/diamond_output/test_XXX_default",
                                                 col_names = F,
                                                 delim="\t",quote = "",comment = "")
CTRLXXX_oxalateproteins_sensitive_raw <- read_delim("~/custom_diamond/diamond_output/test_XXX_sensitive",
                                                   col_names = F,
                                                   delim="\t",quote = "",comment = "")
CTRLXXX_oxalateproteins_moresensitive_raw <- read_delim("~/custom_diamond/diamond_output/test_XXX_moresensitive",
                                                       col_names = F,
                                                       delim="\t",quote = "",comment = "")



#rename columns

new_colnames <- c("query",
                  "subject",
                  "perc_identity",
                  "aln_len",
                  "n_mismatch",
                  "n_gaps",
                  "query_aln_start",
                  "query_aln_end",
                  "subj_aln_start",
                  "subj_aln_end",
                  "evalue",
                  "bitscore")

colnames(CTRLXXX_oxalateproteins_default_raw) <- new_colnames
colnames(CTRXXX_oxalateproteins_sensitive_raw) <- new_colnames
colnames(CTRLXXX_oxalateproteins_moresensitive_raw) <- new_colnames


# sanity check total alignments per query
sanity_check_align <- CTRLXXX_oxalateproteins_default_raw %>% group_by(query) %>% summarize(n=n())
sanity_check_align$n %>% hist(breaks=25) 
# it looks like most queries have either 1 (very unique) or 25 (very common) hits. likely, the nature of the protein database we are using (custom oxalate proteins)
# the max alignment results is 25 but in the graph, we have 50 because we are taking into account the forward and reverse reads as representing/belonging to 1 DNA fragment


# sanity check read pairs ( no longer relevant because we consider fwd and rev reads to be the same)
# sanity_check_pair <- CTRLXXX_oxalateproteins_default_raw %>% group_by(query) %>% slice(1) %>% mutate(read_id = str_extract(query,"HWI.*#") %>% str_replace("#","")) %>% ungroup()
# sanity_check_pair %>% group_by(read_id) %>% summarize(reads = n() %>% factor(levels = c(1,2))) %>% group_by(reads) %>% 
#   summarize(n=n()) %>% ggplot(aes(x=reads,y=n,label=n)) + geom_bar(stat="identity") +geom_label() # about 11% are paired reads
# sanity_check_pair %>% group_by(read_id) %>% mutate(reads = n()) %>% filter(reads==2) %>% mutate(same_query = ifelse(length(unique(subject))==1,TRUE,FALSE)) %>% View()

# add description column to diamond output table
CTRLXXX_oxalateproteins_default_table <- CTRLXXX_oxalateproteins_default_raw %>% 
  inner_join(oxalate_db_map %>% 
               select(subject,description),by="subject") %>% 
  select(query,subject,description,everything())

# sanity check description names
CTRLXXX_oxalateproteins_default_table %>% mutate(nchar = nchar(description)) %>% .$nchar %>% summary()
(CTRLXXX_oxalateproteins_default_table %>% mutate(nchar = nchar(description)) %>% filter(nchar < 827) %>% ggplot(aes(x=nchar)) + geom_histogram(bins=30)) %>% ggplotly
CTRLXXX_oxalateproteins_default_table %>% mutate(nchar = nchar(description)) %>% filter(nchar > 500) %>% View() # wow, there actually exists a description longer than 500 character. I already checked manually,its not an error

# attempt to parse table for abundance
# method: pick first highest bitscore
CTRLXXX_oxalateproteins_default_table$query %>% unique() %>% length() # 44,498 unique query reads

CTRLXXX_oxalateproteins_default_top1 <- CTRLXXX_oxalateproteins_default_table %>% 
  group_by(query) %>% 
  top_n(1,bitscore) %>% 
  slice(1)
CTRLXXX_oxalateproteins_default_top1$query %>% duplicated() %>% which()
CTRLXXX_oxalateproteins_default_top1 %>% View()

CTRLXXX_oxalateproteins_default_top1$description %>% unique() %>% length() # there are 411 unique gene descriptions, too many to plot, we will have to find commonalities between the descriptions
bacteria_extract <- CTRLXXX_oxalateproteins_default_top1 %>% mutate(bacteria_name = str_extract(description,"(\\[)([^\\[]+)(\\]$)"))
bacteria_extract$bacteria_name %>% unique()
bacteria_extract %>% filter(is.na(nchar(bacteria_name))) %>% View()

# To get organism name/lineage from GI or protein name for all proteins, we need to download the xml files (we downloaded fasta files before) and parse it. Might be worth if we cared enough about the organisms
# https://www.biostars.org/p/10959/ or just go to ncbi proteins, search for the protein, and download all proteins in .xml format
nobacteria <- bacteria_extract %>% filter(is.na(nchar(bacteria_name))) 
nobacteria$description %>% unique()

# PUT ON HOLD
# We may have used too many enzymes sequences, making the parsing process more difficult. Some enzymes have rather fixed names, which is great.
# However, other enzymes have complicated names that do not look similar and are very hard to parse. It may be easier to instead use the exact sequences
# that is referenced by curated databases such as KEGG and Metacyc. Currently, I will test this new approach by grabbing the fasta sequences from Metacyc for 
# the gene Oxalate oxidoreductase. See MC8_oxalate_oxidoreductase.R for more

If we redo this anlysis, may want to consider using only uniprot database because the records are better







# # cluster names
# library(RCurl)
# library(stringdist)
# 
# uniquedescriptions <- CTRLXXX_oxalateproteins_default_top1$description %>% unique()
# distancedescriptions <- stringdistmatrix(uniquedescriptions,
#                                          uniquedescriptions,
#                                          method="jw")
# # rownames(distancedescriptions) <- uniquedescriptions
# hc <- hclust(as.dist(distancedescriptions))
# hc
# plot(hc)
# rect.hclust(hc,k=20)
# dfClust <- data.frame(uniquedescriptions, cutree(hc, k=20))
# names(dfClust) <- c('description','cluster')
