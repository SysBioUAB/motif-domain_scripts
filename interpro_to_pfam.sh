#!/bin/bash

# Extract the PFAM records from the InterProScan output. 
if [[ $# != 2 ]] ; then    
    echo 'ERROR: Wrong arguments given. Usage: bash' "$0" '[INPUT_FOLDER_CONTAINING_INTERPROSCAN_OUTPUTS] [OUTPUT_FOLDER]'
    exit
fi
# Change extension if needed ( .fasta or .fasta.tsv)
HOST=$(find $1 -name '*.fasta')

for FILE in $HOST
do
	echo $(basename "$FILE")
	if [[ $(grep "Pfam" $FILE | wc -l) -gt 0 ]]
	then
	  grep "Pfam" $FILE | cut -f 7,8,12,13 > $2$(basename "$FILE")
	fi
done

find $2 -size 0 -delete
