#!/bin/bash
# Make a protein database for the DIAMOND aligner

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

OUTPUT=${OUTPUT}".dmnd"

echo input is $INPUT
echo output is $OUTPUT


# Actual run

if [ "$RUN" = true ]; then

diamond makedb --in "$INPUT" --db "$OUTPUT"

echo Makedb completed!

fi
