#!/bin/bash
# Deduplicate sequences within a protein fasta file based on identical sequence headers

sed -e '/^>/s/$/@@@/' -e 's/^>/###/' custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_linear.fasta  | tr -d '\n' | sed -e 's/###/\'$'\n/g'  -e 's/@@@/\'$'\t/g' | sort -u -t$'\t' -f -k1,1 | sed -e 's/^/>/' | tr "\t" "\n" > ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique.fasta

# remove empty first line
sed -e '1!b' -e '/^>gi/!d' ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique.fasta > ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique2.fasta

# add newline to separate sequences. this fixed makedb error for diamond
sed -e '/^>gi/b' -e 's/$/\'$'\n/' ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique2.fasta > ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique.fasta

# remove temp file
rm ~/custom_diamond/db/oxalate_enzymes_ncbi_alldb/DB_oxalate_proteins_unique2.fasta

echo Protein sequences deduplicated!