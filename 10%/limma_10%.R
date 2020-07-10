if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

options(repos = c(
  "https://iwww.broadinstitute.org/~datasci/R-packages",
  "https://cran.cnr.berkeley.edu"))
#install.packages('cdsr')

library(taigr)

#gene.effect <- load.from.taiga(data.name='avana-public-tentative-19q4-c2df', data.version=4, data.file='gene_effect')
gene.effect <-  load.from.taiga(data.name='demeter2-achilles-5386', data.version=13, data.file='gene_effect')
#gene.effect <- load.from.taiga(data.name='demeter2-drive-0591', data.version=12, data.file='gene_effect')

aneuploidy.data <- load.from.taiga(data.name='aneuploidy-data-d0b9', data.version=19, data.file='using_ten_percent_NEW')

high_aneuploidy = aneuploidy.data[,'many_arm_events'] == 'True'
low_aneuploidy = aneuploidy.data[,'many_arm_events'] == 'False'

names_high = aneuploidy.data[high_aneuploidy,'CCLE_ID']
names_low = aneuploidy.data[low_aneuploidy,'CCLE_ID']


both = intersect(rownames(gene.effect),aneuploidy.data[,'CCLE_ID'])
#both = intersect(rownames(gene.expression),aneuploidy.data[,'DepMap_ID'])

gene.effect_2 = gene.effect[both,]

in_group = row.names(gene.effect_2) %in% names_high

# if including as a covariate
#sample.info <- load.from.taiga(data.name='public-19q4-93d9', data.version=21, data.file='sample_info')

source("running_limma.R")

lim_res <- run_lm_stats_limma(gene.effect_2, in_group)#,covar = sample.info[rownames(gene.expression_2),'lineage',drop=F])
write.table(lim_res,file='lim_res_rnai_achilles_10_percent_NEW.csv',sep = ",",row.names = FALSE)



