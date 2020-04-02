#!/bin/bash
# Deduplicate sequences within a protein fasta file based on identical sequence headers


# Options and variables

RUN=false

while getopts i:d:o:r option
do
        case "${option}"
        in
                i) INPUT=${OPTARG};;
                d) INPUT_DIR=${OPTARG};;
                o) OUTPUT=${OPTARG};;
		            r) RUN=true;;
        esac
done

INPUT_DIR=${INPUT_DIR}"/"
TEMP1=${INPUT_DIR}"temp_unique1.fasta"
TEMP2=${INPUT_DIR}"temp_unique2.fasta"

echo input is $INPUT
echo input_dir is $INPUT_DIR
echo temp1 is $TEMP1
echo temp1 is $TEMP2
echo output is $OUTPUT 


if [ "$RUN" = true ]; then

sed -e '/^>/s/$/@@@/' -e 's/^>/###/' "$INPUT"  | tr -d '\n' | sed -e 's/###/\'$'\n/g'  -e 's/@@@/\'$'\t/g' | sort -u -t$'\t' -f -k1,1 | sed -e 's/^/>/' | tr "\t" "\n" > "$TEMP1"

# remove empty first line
sed -e '1!b' -e '/^>gi/!d' "$TEMP1" > "$TEMP2"

# add newline to separate sequences. this fixed makedb error for diamond
sed -e '/^>gi/b' -e 's/$/\'$'\n/' "$TEMP2" > "$OUTPUT"

# remove temp file
rm "$TEMP1" 
rm "$TEMP2"

echo Protein sequences deduplicated!

fi






