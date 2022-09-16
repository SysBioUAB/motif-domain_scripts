import os.path
import glob
import re
import os
import pandas as pd
import subprocess

# List all the human proteins with the windowed accessibility scores from processedalphafold github.
host_files = [os.path.basename(x).split(".")[0] for x in glob.glob("./host_windowed_scores/*.txt")]

# Dataframe obtained in the clean_phisto.py with all the host-pathogen PPI
df = pd.read_csv('host-pathogen_interaction.csv')

# Function that takes the accessibility scores, if the score is higher than 0.55 (best option according to the authors)
# and if there are 5 or more contiguous residues, the region is stored as IDR. The positions and the sequence.
def host_IDR(files):
    for file in files:
        with open("./host_windowed_scores/"+file+'.txt','rt') as f:
            with open("./host_fasta/"+file+".fasta",'rt') as f_in:
                sequence = f_in.readlines()[1]
                accessibility_score = f.readline().split('	')[1].split(',')
                positions = []
                for i in range(len(accessibility_score)):
                    if float(accessibility_score[i]) >= 0.55:
                        positions.append(i)
                    else:
                        if len(positions) > 5:
                            with open("./host_IDR/"+file+'.txt','at') as f_out:
                                #print(str(positions[0])+'\t'+str(positions[-1])+'\t'+sequence[positions[0]:positions[-1]])
                                f_out.write(str(positions[0])+'\t'+str(positions[-1])+'\t'+sequence[positions[0]:positions[-1]]+'\n')

                        positions = []

# This function requires the ELM motif database downloaded and with only the columns corresponding to the motif ID and the regex.
# For each protein that contain IDR and for each IDR, search if there are motifs and stores it in a new file.
def search_motifs(files):
    regex = []
    ELM = []

    with open('ELM_motifs.tsv', 'rt') as f:
        for line in f.readlines():
            if line[0] != '#' and line[1] == 'E':
                ELM.append(line.split("\t")[0].strip().split('"')[1])
                regex.append(line.split("\t")[4].strip().split('"')[1])

    # This change in the regex is necessary for the findall function.
    # The problem  is that if the regex that re.findall tries to match captures groups
    # (i.e. the portions of the regex that are enclosed in parentheses), then it is the groups that are returned, rather than the matched string.
    regex = [i.replace('(', '(?:') for i in regex]

    for file in files:
        try:
            with open("./host_IDR/"+file+'.txt','rt') as f:
                for pos_IDR in f.readlines():
                    for i in range(len(regex)):

                        # Instead of using the re.search and re.groups, the re.compile, re.findall and re.finditer are used
                        r = re.compile(regex[i])
                        matches = re.findall(regex[i],pos_IDR.split('\t')[-1])
                        if len(matches) > 0:
                            for j in range(len(matches)):
                                if isinstance(matches[j],tuple):
                                    for match in matches[j]:
                                        if len(match) >= 5:
                                            start_end = [[m.start(), m.end()] for m in r.finditer(pos_IDR.split('\t')[-1])][j]
                                            with open("./host_motif/" + file + '.txt', 'at') as f:
                                                f.write("%s\t%s\t%s\t%s\n" % (
                                                str(int(pos_IDR.split('\t')[0]) + int(start_end[0])),
                                                str(int(pos_IDR.split('\t')[0]) + int(start_end[1])), ELM[i],
                                                str(matches[j])))

                                elif len(matches[j]) >= 5 and isinstance(matches[j],str):
                                    start_end = [[m.start(), m.end()] for m in r.finditer(pos_IDR.split('\t')[-1])][j]
                                    with open("./host_motif/"+file+'.txt','at') as f:
                                        f.write("%s\t%s\t%s\t%s\n" % (str(int(pos_IDR.split('\t')[0])+int(start_end[0])),str(int(pos_IDR.split('\t')[0])+int(start_end[1])),ELM[i],str(matches[j])))
        except:
            with open("./host_motif/" + file + '.txt', 'at') as f:
                f.write("")

# Function to generate all the motif-domain (host-pathogen) combinations.
def motif_domain_combinations(df):
    for i in df.itertuples():
        print(i._3 + '-' + i._2)
        try:
            with open('host_motif/' + i._3 + '.txt', 'rt') as f_host:
                motifs = f_host.readlines()
                with open('Pfam_domains/' + i._2 + '.fasta.tsv', 'rt') as f_pathogen:
                    domains = f_pathogen.readlines()
                    with open('motif_domain_combinations/' + i._3 + '-' + i._2 + '.txt', 'wt') as f:
                        f.write('Host_uniprot: ' + i._3 + '\nPathogen_uniprot: ' + i._2 + '\nPathogen: ' + \
                                i.Pathogen + '\nHost: Homo sapiens' + '\nHost_sequence: ' + i.Sequence_host \
                                + '\nPathogen_sequence: ' + i.Sequence_pathogen + '\n')
                        for host_line in motifs:
                            for pathogen_line in domains:
                                if host_line.strip() != '':
                                    f.write(host_line + pathogen_line + '\n\n')
        except:
            continue
            #print('Files not found.')

# host_IDR(host_files)
#search_motifs(host_files)
motif_domain_combinations(df)


