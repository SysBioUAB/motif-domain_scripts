library(rstatix)
library(coin)

# Load all the tables (better using the import dataset text to select the separartors and look at the data)

# Initial dataframe with headers 
results <- data.frame(domain = 'Pfam ID', wilcoxon = 'Wilcoxon p-value', effect_size= 'Effect Size')

# For each domain compute a one-sample Wilcoxon test using as a data vector the frequencies of the random bootstrapping and as a mu
# the frequency of the observed domain.  
for (i in seq(1,nrow(domain.domain_random_freq_prueba))){
  
  # The name of the domain as a character
  domain <- as.character(domain.domain_random_freq_prueba[i,1])
  
  # The observed frequency as a number. The decimal was separated with comma, now with a dot. This is going to be the mu of the wilcoxon test.
  freq <- as.numeric(sub(",",".",domain.domain_observed_freq_prueba[domain.domain_observed_freq_prueba$V2 == domain,3],fixed = TRUE))
  
  # From the random frequencies, take the row where the first column is the domain (dataframe containing 1 row with 1000 columns with frequencies)
  distri <- data.frame(domain.domain_random_freq_prueba[domain.domain_random_freq_prueba$V1 == domain,c(-1)])
  
  # Transpose matrix so it has 1 column and 1000 rows, change the column name and transform it into a dataframe, so the
  # wilcox_effsize function can work properly. The %>% is a pipe, domain_freq is the name of the column, the 1 is used
  # because it is a One-sample wilcoxon test (for a two-sample it would look like domain_freq ~ another_column_name), and the mu
  # forms the null hypothesis.
  transpose <- t(distri)
  colnames(transpose) <- c("domain_freq")
  domain_effect_size <- as.data.frame(transpose) %>% wilcox_effsize(domain_freq ~ 1,mu=freq)
  
  # Computes the wilcoxon test.
  test <- wilcox.test(as.numeric(distri), mu=freq)
  
  # It creates a dataframe with the same shape as results, so the rows can be bound. By default it creates 1 column and 3 rows 
  # (domain, wilcoxon, effect_size) and with the transposition, it is 1 row with 3 columns, as in the results dataframe.
  # The column name of the effect_size is changed so it is the same as in results (Using the domain_effect_size$effsize the column name)
  # is modified somehow, so we correct it. 
  newline <- data.frame(t(c(domain= domain,wilcoxon = test$p.value, effect_size = domain_effect_size$effsize)))
  colnames(newline) <- c('domain','wilcoxon','effect_size')
  results <- rbind(results,newline)
}

# At the end of the loop, all the results are written in a csv file.
write.table(data.frame(results),file = '~/Jordi/Host_pathogen_project_copy/Host_pathogen/outputs/domain-domain_wilcoxon_test.txt', row.names = FALSE, append = FALSE, col.names = FALSE, sep = ", ")

library('ggpubr')

host_motif_enrichment <- read.table('~/Jordi/Host_pathogen_project_copy/Host_pathogen/outputs/host_motif_enrichment.txt')
host_domain_enrichment <- read.table('~/Jordi/Host_pathogen_project_copy/Host_pathogen/outputs/host_domain_enrichment_01.txt') 
pathogen_domain_enrichment <- read.table('~/Jordi/Host_pathogen_project_copy/Host_pathogen/outputs/pathogen_domain_enrichment_05.txt')
domain_domain_enrichment <- read.table('~/Jordi/Host_pathogen_project_copy/Host_pathogen/outputs/domain-domain_enrichment_nuevo.txt')
motif_domain_enrichment <- read.table('~/Jordi/Host_pathogen_project_copy/Host_pathogen/outputs/motif-domain_enrichment_01_nuevo.txt')
  
df <- motif_domain_enrichment
pos_ordered <- c(order(df$V2))
pos <- tail(pos_ordered,20)

df_top10 <- df[pos,]

df_top10$grp <- factor(ifelse(df_top10$V2 < 1.0, "red", "black"), 
                 levels = c("red", "black"))
colnames(df_top10) <- c('domain','frequency','color')
barplot <- ggbarplot(df_top10, x = "domain", y = "frequency",
          title="Top 20 enriched motif-domain associations",
          fill=df_top10$color,
          color = "white",   
          label= F,
          lab.col="black",
          lab.nb.digits = 4,
          lab.hjust = 0.5,
          sort.by.groups = FALSE,     
          lab.vjust = -1,
          x.text.angle = 90,    
          xlab = "Motif-domain",
          ylab = "Fobs/Fexp",
          rotate = TRUE,
          ggtheme = theme_minimal()
) + theme(plot.title = element_text(hjust=0.4))
                                                                     
ggpar(barplot,font.main = c(16,"bold"), font.xtickslab = c(16,"bold"), font.ytickslab = c(16,"bold"),font.x = c(16,"bold"),font.y = c(16,"bold")) +coord_cartesian(ylim = c(0,15)) 




