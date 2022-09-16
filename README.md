# Analysis of host-bacteria protein interactions reveals conserved domains and motifs that mediate fundamental infection pathways

These scripts allow for the preprocessing steps of the PHISTO database and the generation of the motif/domain combinations.

clean_phisto.py: Clean the PHISTO database and generate the FASTA files required to run InterProScan. Generates all the domain combinations given the InterProScan output. Multiple functions are defined in this script: 
- clean_phisto() takes as an input the PHISTO database (phi_data.csv) and outputs the UniProt IDs of all the host/pathogen proteins. 
- sequences_and_fasta() generates the FASTA files required for InterProScan to predict the domains.
- domain_combinations() generates all the domain combinations between the host and pathogen proteins.

interpro_to_pfam.sh: Extract only the PFAM records from the Interproscan output. Usage: bash interpro_to_pfam.sh [INPUT_FOLDER_CONTAINING_INTERPROSCAN_OUTPUTS] [OUTPUT_FOLDER]'. The host_fasta and pathogen_fasta folders are used as input and contain the outputs from InterProScan (.fasta.tsv files). InterProScan v5.56 was used (bash interproscan.sh -i fasta_file -f tsv).

host_motifs.py: Compute the IDR regions of the human proteins and find the ELM motifs in these regions. Generates the MD and DD combinations.

compute_freq.sh: Compute the relative frequencies.

LICENSE

![image](https://user-images.githubusercontent.com/78474998/190620662-51e972db-df9d-42cf-a215-758c58e9e5f3.png)

Creative Commons License
This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
