#!/bin/bash
# Make a protein database for the DIAMOND aligner

diamond makedb --in ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique.fasta --db ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique

echo Makedb completed!