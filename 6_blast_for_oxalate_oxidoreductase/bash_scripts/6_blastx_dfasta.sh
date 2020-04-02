#!/bin/bash
# bash command to run diamond's blastx (translated DNA query aligned to protein db)


# Options and variables

RUN=false
SETTING="default"
SETTING_CASE=""

while getopts i:d:o:s:r option
do
        case "${option}"
        in
                i) INPUT=${OPTARG};;
                d) DB=${OPTARG};;
                o) OUTPUT=${OPTARG};;
                s) SETTING=${OPTARG};;
                r) RUN=true;;
        esac
done

case $SETTING
        in
                "default" ) SETTING_CASE="";;
                "sensitive" ) SETTING_CASE="sensitive";;
                "more-sensitive" ) SETTING_CASE="more-sensitive";;
esac

# uncomment for details of run
# echo input is: $INPUT
# echo database is: $DB
# echo output is: $OUTPUT
# echo setting is: $SETTING
# echo command is: diamond blastx --db "$DB" -q "$INPUT" ${SETTING_CASE:+ "--$SETTING_CASE"} -o "$OUTPUT"


# Actual run

if [ "$RUN" = true ]; then

diamond blastx --db "$DB" -q "$INPUT" ${SETTING_CASE:+ "--$SETTING_CASE"} -o "$OUTPUT"

echo Blastx complete!

fi
