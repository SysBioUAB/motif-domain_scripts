#!/bin/bash

# For bootstrap_host
if [[ $# == 1 ]]
then  
    # Counting the number of unique domains
    DOMAIN_COUNTS=$(cat $1 | cut -f3 | sort | uniq -c)
fi

# For bootstrap_pathogen
if [[ $# == 3 ]]
then  
    # Counting the number of unique domains
    DOMAIN_COUNTS=$(cat $1 $2 $3 | cut -f3 | sort | uniq -c)
fi

# Taking the total number of domains
    TOTAL=$(echo "$DOMAIN_COUNTS" | awk '{ sum += $1 } END { print sum }')

    # Compute the frequency for all the domains 
    echo "$DOMAIN_COUNTS" | awk -v sum="$TOTAL" '{ print $0, $1/sum }' 

