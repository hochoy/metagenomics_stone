library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(R.utils)

# PRE. first step is to search for all the proteins of interest in NCBI and download all (I did this manually but we could use NCBI's Entrez API)
# Since we are looking for similar proteins, it's possible that different queries yield similar proteins.
# E.g. A search for Pyruvate oxidoreductase vs pyruvate decarboxylase may return similar proteins.
# In that case, we need to remove duplicate sequences that share the same FASTA sequence header. 
# Additionally, we could also cluster "similar" sequences using CD-HIT or ClustalW to further reduce the number of proteins (not done).


# 1.Concatenate the protein.fasta files
# cat ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/source/*.fasta > ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins.fasta
# system("~/aNALYSIS_PIPELINE/git_world/gut_metagenomics/metagenomics/MC8_oxalateenzymes/old/bash_scripts/1_cat_pfasta.sh")


# 2. linearize protein fasta file (protein sequences are currently wrapped: meaning there's a newline separating each protein seq line)
# sed -e 's/\(^>.*$\)/###\1###/' ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins.fasta | tr -d "\r" | tr -d "\n" | sed -e 's/$/###/' | sed -e 's/###/\'$'\n/g' | sed -e '/^$/d' > ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_linear.fasta
# system("~/aNALYSIS_PIPELINE/git_world/gut_metagenomics/metagenomics/MC8_oxalateenzymes/old/bash_scripts/2_linearize_pfasta.sh")


# 3. deduplicate protein sequences (two or more sequences that have identical headers)
# sed -e '/^>/s/$/@@@/' -e 's/^>/###/' custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_linear.fasta  | tr -d '\n' | sed -e 's/###/\'$'\n/g'  -e 's/@@@/\'$'\t/g' | sort -u -t$'\t' -f -k1,1 | sed -e 's/^/>/' | tr "\t" "\n" > ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique.fasta
# sed -e '1!b' -e '/^>gi/!d' ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique.fasta > ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique2.fasta
# sed -e '/^>gi/b' -e 's/$/\'$'\n/' ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique2.fasta > ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique.fasta
# rm ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique2.fasta
# system("~/aNALYSIS_PIPELINE/git_world/gut_metagenomics/metagenomics/MC8_oxalateenzymes/old/bash_scripts/3_unique_pfasta.sh")


# 4. grep protein fasta headers
# grep ">" ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique.fasta > ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_header.fasta
# system("~/aNALYSIS_PIPELINE/git_world/gut_metagenomics/metagenomics/MC8_oxalateenzymes/old/bash_scripts/4_header_pfasta.sh")


# 5. makedb
# diamond makedb --in ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique.fasta --db ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique
# system("~/aNALYSIS_PIPELINE/git_world/gut_metagenomics/metagenomics/MC8_oxalateenzymes/old/bash_scripts/5_makedb_pfasta.sh")


# 6. run diamond
# diamond blastx --db ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique -q ~/rawread_qc/rawreads_qced/humann2_input/CTRLXXX_merged.fastq.gz -o ~/custom_diamond/diamond_output/test_A6_default
# diamond blastx --db ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique -q ~/rawread_qc/rawreads_qced/humann2_input/CTRLXXX_merged.fastq.gz --sensitive -o ~/custom_diamond/diamond_output/test_A6_sensitive
# diamond blastx --db ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique -q ~/rawread_qc/rawreads_qced/humann2_input/CTRLXXX_merged.fastq.gz --more-sensitive -o ~/custom_diamond/diamond_output/test_A6_moresensitive
# system("~/aNALYSIS_PIPELINE/git_world/gut_metagenomics/metagenomics/MC8_oxalateenzymes/old/bash_scripts/6_blastx_dfasta.sh")

