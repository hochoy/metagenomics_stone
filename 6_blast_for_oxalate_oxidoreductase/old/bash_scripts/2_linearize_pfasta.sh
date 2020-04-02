#!/bin/bash
# Linearize a protein fasta file

sed -e 's/\(^>.*$\)/###\1###/' ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins.fasta | tr -d "\r" | tr -d "\n" | sed -e 's/$/###/' | sed -e 's/###/\'$'\n/g' | sed -e '/^$/d' > ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_linear.fasta



# Original comments: A useful step is to linearize your sequences (i.e. remove the sequence wrapping). This is not a perfect solution, as I suspect that a few steps could be avoided, but it works quite fast, even for thousands of sequences.

echo Fasta file linearized!