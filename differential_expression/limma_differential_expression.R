if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

options(repos = c(
  "https://iwww.broadinstitute.org/~datasci/R-packages",
  "https://cran.cnr.berkeley.edu"))
#install.packages('cdsr')

library(taigr)


gene.expression <- load.from.taiga(data.name='public-19q4-93d9', data.version=21, data.file='CCLE_expression')

aneuploidy.data <- load.from.taiga(data.name='aneuploidy-data-d0b9', data.version=18, data.file='aneuploidy_data_NEW_COMPACT')

high_aneuploidy = aneuploidy.data[,'many_arm_events'] == TRUE
low_aneuploidy = aneuploidy.data[,'many_arm_events'] == FALSE

names_high = aneuploidy.data[high_aneuploidy,'DepMap_ID']
names_low = aneuploidy.data[low_aneuploidy,'DepMap_ID']


both = intersect(rownames(gene.expression),aneuploidy.data[,'DepMap_ID'])

gene.expression_2 = gene.expression[both,]

in_group = row.names(gene.expression_2) %in% names_high


sample.info <- load.from.taiga(data.name='public-19q4-93d9', data.version=21, data.file='sample_info')
rownames(sample.info) <- sample.info$DepMap_ID #if using depmap IDs

source("running_limma.R")

## uncomment if including lineage as a covariate
lim_res <- run_lm_stats_limma(gene.expression_2, in_group)#,covar = sample.info[rownames(gene.expression_2),'lineage',drop=F])


write.table(lim_res,file='differential_expression_NEW.csv',sep = ",",row.names = FALSE)



