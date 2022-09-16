import sys
# Domain-domain combinations for the bootstrapping domain-domain_bootstrap.sh.
iteration = sys.argv[1]

with open('DOMAIN_COMBINATION'+iteration+'.txt','rt') as f:
    for line in f:
        host = line.strip().split()[1]
        pathogen = line.strip().split()[2]
        with open(host,'rt') as f_host, open(pathogen, 'rt') as f_path:
            host_domains = f_host.readlines()
            pathogen_domains = f_path.readlines() 
            for host_domain in host_domains:
            	#uncomment for motif-domain
            	#if len(host_domain) > 1:
            		for pathogen_domain in pathogen_domains:
            			print(host_domain.split()[2]+'\t'+pathogen_domain.split()[2])
