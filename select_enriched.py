from matplotlib import pyplot as plt
import subprocess
import statistics
import os

host_pathogen = 'domain'
domain_or_motif = 'domain'
cutoff = 5

interactions_with_obs_freq = {}
interactions_with_rand_freq = {}

# From the results of the wilcoxon test and effect size, take domains where p-value > 0.05 and effect size > 0.5
cat = subprocess.Popen('cat outputs/' + host_pathogen+'-'+domain_or_motif+'_wilcoxon_test.txt', stdout=subprocess.PIPE, shell=True)
tr = subprocess.Popen('tr \'\",\' \'  \'', stdin=cat.stdout, stdout=subprocess.PIPE, shell=True)
awk = subprocess.Popen('LC_NUMERIC=C awk \'{if ($3 + 0 > 0.5 && $2 + 0 < 0.05) {print $1}}\' ',stdout=subprocess.PIPE,stdin=tr.stdout,shell=True)

domains = awk.communicate()[0].decode().split('\n')[1:]

# Take only the enriched motif/domains from the host/pathogen (values above 1)
motif_enrichment = subprocess.Popen('LC_NUMERIC=C awk \'{if ($2 + 0 > 1) {print $1}}\' outputs/host_domain_enrichment.txt',stdout=subprocess.PIPE,shell=True)
motif_enrichment_only_enriched = motif_enrichment.communicate()[0].decode().split('\n')[1:]

domain_enrichment = subprocess.Popen('LC_NUMERIC=C awk \'{if ($2 + 0 > 1) {print $1}}\' outputs/pathogen_domain_enrichment.txt',stdout=subprocess.PIPE,shell=True)
domain_enrichment_only_enriched = domain_enrichment.communicate()[0].decode().split('\n')[1:]


with open('outputs/domain-domain_observed_freq.txt','rt') as f_obs:
	for line in f_obs:
		record_obs = line.split()
		if int(record_obs[0]) >= cutoff:
			interactions_with_obs_freq[record_obs[-2]] = record_obs[-1].replace(',','.')


with open('outputs/domain-domain_random_freq.txt', 'rt') as f_rand:
	for line in f_rand:
		record_rand = line.split(',')
		all_rand_freq = record_rand[1:]
		rand_freq = statistics.fmean([float(i) for i in all_rand_freq])
		interactions_with_rand_freq[record_rand[0]] = rand_freq

# Compute the fobs/frand for the selected domains.
for i in domains:
	try:
		if i in interactions_with_obs_freq:
			motif = i.split('-')[0]
			domain = i.split('-')[1]
			n=0

			if motif in motif_enrichment_only_enriched and domain in domain_enrichment_only_enriched:

				prot1 = subprocess.Popen("grep -nrw " + motif + " Pfam_domains_combinations/", stdout=subprocess.PIPE,
											shell=True)
				prot2 = subprocess.Popen("grep -nrw " + domain + " Pfam_domains_combinations/", stdout=subprocess.PIPE,
											shell=True)

				#num_prot = subprocess.Popen("grep -nrw " + i  + " Pfam_domains_pathogen/", stdout=subprocess.PIPE, shell=True)
				protein_with_motif = prot1.communicate()[0].decode().split("\n")
				protein_with_domain = prot2.communicate()[0].decode().split("\n")

				unique_motif = []
				unique_domain = []

				for j in protein_with_motif:
					if j != '' and j.split('/')[1].split('.')[0] not in unique_motif:
						unique_motif.append(j.split('/')[1].split('.')[0])
				for j in protein_with_domain:
					if j != '' and j.split('/')[1].split('.')[0] not in unique_domain:
						unique_domain.append(j.split('/')[1].split('.')[0])
				for j in set(unique_domain):
					if j in unique_motif:
						print(j)
						n+=1

				if i in interactions_with_rand_freq and n >= 1:

					obs_freq = float(interactions_with_obs_freq[i])
					rand_freq = float(interactions_with_rand_freq[i])


					with open('outputs/domain_domain_enrichment.txt','at') as f:
						f.write(i+'\t'+str(obs_freq/rand_freq)+'\n')

	except:
		continue
