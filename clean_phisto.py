import pandas as pd
import subprocess

# This function clean the phisto database and create files with uniprot ID for host and pathogen to use the Uniprot mapper.
def clean_phisto():
    # Loading data and cleaning
    df = pd.read_csv('phi_data.csv')
    df_clean = df[['Pathogen','Uniprot ID','Uniprot ID.1']]

    # Uniprot ids in lists
    pathogen_uniprot = df_clean['Uniprot ID'].tolist()
    host_uniprot = df_clean['Uniprot ID.1'].tolist()

    # Create files with uniprot IDs to obtain sequences in uniprot mapper.
    with open('pathogen_uniprot.txt','wt') as f_pat:
        with open('host_uniprot.txt', 'wt') as f_host:
            for i in range(len(pathogen_uniprot)):
                f_pat.write(pathogen_uniprot[i]+'\n')
                f_host.write(host_uniprot[i]+'\n')

    return df_clean
# Once the first function is executed, the uniprot ids are mapped in uniprot mapper and the excel is downloaded and located at the working directory.
# This function requires the clean df from the previous function and the excels from uniprot mapper for host and pathogen (uniprot_host.xlsx and uniprot_pathogen.xlsx )
# These df are cleaned and all the FASTA files needed for the interproscan are generated.
def sequences_and_fasta(df_clean):
    # Read excels with sequences and clean data
    df_host = pd.read_excel('uniprot_host.xlsx')
    df_pathogen = pd.read_excel('uniprot_pathogen.xlsx')
    df_host_clean = df_host[['Entry','Sequence']]
    df_pathogen_clean = df_pathogen[['Entry','Sequence']]
    df_host_clean.columns = ['Uniprot ID.1','Sequence_host']
    df_pathogen_clean.columns = ['Uniprot ID','Sequence_pathogen']

    # Merge main dataframe with sequences
    merge_host_seq = df_clean.merge(df_host_clean, on='Uniprot ID.1',how='left')
    merge_pathogen_seq = merge_host_seq.merge(df_pathogen_clean,on='Uniprot ID', how='left')
    df_final = merge_pathogen_seq
    df_final.to_csv('host-pathogen_interaction.csv', index=False)

    # Generate all the fasta files

    for i in df_pathogen_clean.itertuples():
        with open('pathogen_fasta/'+i._1+'.fasta', 'wt') as f:
            f.write('>'+i._1+'\n'+i.Sequence_pathogen)

    for i in df_host_clean.itertuples():
        with open('host_fasta/'+i._1+'.fasta', 'wt') as f:
            f.write('>'+i._1+'\n'+i.Sequence_host)
    return df_final

# Function to generate all the domain combinations between the host and pathogen proteins
def domain_combinations(df):
    for i in df.itertuples():
        print(i._3 + '-' + i._2)
        try:
            with open('Pfam_domains/' + i._3 + '.fasta.tsv', 'rt') as f_host:
                host = f_host.readlines()
                with open('Pfam_domains/' + i._2 + '.fasta.tsv', 'rt') as f_pathogen:
                    pathogen = f_pathogen.readlines()
                    with open('Pfam_domains_combinations/' + i._3 + '-' + i._2 + '.txt', 'wt') as f:
                        f.write('Host_uniprot: ' + i._3 + '\nPathogen_uniprot: ' + i._2 + '\nPathogen: ' + \
                                i.Pathogen + '\nHost: Homo sapiens' + '\nHost_sequence: ' + i.Sequence_host \
                                + '\nPathogen_sequence: ' + i.Sequence_pathogen + '\n')
                        for host_line in host:
                            for pathogen_line in pathogen:
                                f.write(host_line + pathogen_line + '\n\n')

        except:
            continue


# Step 2
df_clean = clean_phisto()
# Step 4
df_final = sequences_and_fasta(df_clean)
# Step 6
# subprocess.run(["bash", "extract_pfam_from_tsv.sh"])
# Step 7
print(domain_combinations(df_final))
