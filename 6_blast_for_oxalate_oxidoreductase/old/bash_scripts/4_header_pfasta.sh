#!/bin/bash
# Grep sequence headers from a protein fasta file

grep ">" ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique.fasta > ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_header.fasta

echo Fasta headers grepped!