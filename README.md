# domain_scripts

These scripts allows for the preprocessing steps of the PHISTO database and the generation of the motif/domain combinations.

clean_phisto.py: Clean the PHISTO database and generate the FASTA files required to run InterProScan. Also generates all the domain combinations given the InterProScan output.

host_motifs.py: Compute the IDR regions of the human proteins and fetch the ELM motifs in these regions. Also generates the MD and DD combinations.

compute_freq.sh: Compute the relative frequencies of the domains.

LICENSE

Creative Commons License
This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
