#!/bin/bash
# Grep sequence headers from a protein fasta file


# Options and variables

RUN=false

while getopts i:o:r option
do
        case "${option}"
        in
                i) INPUT=${OPTARG};;
                o) OUTPUT=${OPTARG};;
                r) RUN=true;;
        esac
done

echo input is $INPUT
echo output is $OUTPUT

# Actual run

if [ "$RUN" = true ]; then

grep ">" "$INPUT" > "$OUTPUT"

echo Fasta headers grepped!

fi
