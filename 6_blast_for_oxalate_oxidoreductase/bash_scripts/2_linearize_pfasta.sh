#!/bin/bash
# Linearize a protein fasta file

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

sed -e 's/\(^>.*$\)/###\1###/' "$INPUT" | tr -d "\r" | tr -d "\n" | sed -e 's/$/###/' | sed -e 's/###/\'$'\n/g' | sed -e '/^$/d' > "$OUTPUT"

fi

# Original comments: A useful step is to linearize your sequences (i.e. remove the sequence wrapping). This is not a perfect solution, as I suspect that a few steps could be avoided, but it works quite fast, even for thousands of sequences.

echo Fasta file linearized!


