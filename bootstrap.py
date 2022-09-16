from joblib import Parallel, delayed
from matplotlib import pyplot as plt
import scipy.stats
import subprocess
import os

# Number of bootstrap iterations and dictionary containing domain frequencies.
pfam_paths = []
steps = 1000
domain_freq = {}
host_or_pathogen = 'host'
domain_or_motif = 'motif'

if host_or_pathogen + '_random_'+domain_or_motif+'_freq.txt' in os.listdir('./outputs/'):
    os.remove('outputs/'+host_or_pathogen+'_random_'+domain_or_motif+'_freq.txt')

# Function to execute the bootstrap script in parallel.
# host: '../Pfam/human_Pfam/'
# pathogen : '../Pfam/yersinia_Pfam/' ,'../Pfam/bacillus_Pfam/ ','../Pfam/francisella_Pfam/'

def bootstrap():
    result = subprocess.run(
        ['bash', 'bootstrap_' + host_or_pathogen + '.sh', '../Pfam/human_motif/'], stdout=subprocess.PIPE)
    return result.stdout.decode().split('\n')

# Calling the previous function in parallel, adjust the number of jobs if needed (maybe 6 is better) and put the pfam_path.
par = Parallel(n_jobs=4)(delayed(bootstrap)() for i in range(steps))


# Stores in the dictionary the names of the domains in the keys and the frequencies in the values as a list.
for iteration in par:
    for domain in iteration:
        pfam_freq = domain.split(' ')[-2:]
        if pfam_freq != ['']:
            if pfam_freq[0] not in domain_freq:
                domain_freq[pfam_freq[0]] = [float(pfam_freq[1].replace(',','.'))]
            else:
                domain_freq[pfam_freq[0]] += [float(pfam_freq[1].replace(',','.'))]


# We are only interested in the domains that we observed from the phisto database, to the rest of domains obtained
# at random from the proteomes is not useful. A shapiro test is computed to know if the frequency distributions are
# normal. A file with the name of the domain and the list of frequencies is created.
print('yes')
with open('outputs/'+host_or_pathogen+'_observed_'+domain_or_motif+'_freq.txt','rt') as f_obs:
    lines = f_obs.readlines()
    for i in domain_freq.items():
        for j in lines:
            if i[0] in j and i[0] != '':
                freq_with_zeros = i[1] + [0 for i in range(steps-len(i[1]))]
                #shapiro = scipy.stats.shapiro(freq_with_zeros)
                with open('outputs/'+host_or_pathogen+'_random_'+domain_or_motif+'_freq.txt','at') as f_rand:
                    #f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (i[0],str(freq_with_zeros)[1:len(str(freq_with_zeros))-1], shapiro.statistic,shapiro.pvalue,w,p))
                    f_rand.write('%s,%s\n' % (i[0], str(freq_with_zeros)[1:len(str(freq_with_zeros)) - 1]))

