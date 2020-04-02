#!/bin/bash
# bash command to run diamond's blastx (translated DNA query aligned to protein db)

diamond blastx --db ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique -q ~/rawread_qc/rawreads_qced/humann2_input/CTRLXXX_merged.fastq.gz -o ~/custom_diamond/diamond_output/test_A6_default

diamond blastx --db ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique -q ~/rawread_qc/rawreads_qced/humann2_input/CTRLXXX_merged.fastq.gz --sensitive -o ~/custom_diamond/diamond_output/test_A6_sensitive

diamond blastx --db ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique -q ~/rawread_qc/rawreads_qced/humann2_input/CTRLXXX_merged.fastq.gz --more-sensitive -o ~/custom_diamond/diamond_output/test_A6_moresensitive

echo Blastx complete!