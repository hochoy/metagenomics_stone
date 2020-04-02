#!/bin/bash
# Combine all the oxalate protein fastas

cat ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/source/*.fasta > ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins.fasta

echo Catenated oxalate proteins!