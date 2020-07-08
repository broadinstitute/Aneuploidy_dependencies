if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

options(repos = c(
  "https://iwww.broadinstitute.org/~datasci/R-packages",
  "https://cran.cnr.berkeley.edu"))
install.packages('cdsr')

library(taigr)

# get prism secondary data
secondary.matrix <- load.from.taiga(data.name='aneuploidy-data-d0b9', data.version=5, data.file='secondary_matrix')
aneuploidy.data <- load.from.taiga(data.name='aneuploidy-data-d0b9', data.version=18, data.file='aneuploidy_data_NEW_COMPACT')


high_aneuploidy = aneuploidy.data[,'many_arm_events'] == TRUE
low_aneuploidy = aneuploidy.data[,'many_arm_events'] == FALSE

both = intersect(rownames(secondary.matrix),aneuploidy.data[,'DepMap_ID'])
sub = secondary.matrix[both,]

names_high = aneuploidy.data[high_aneuploidy,'DepMap_ID'] # for crispr
in_group = row.names(sub) %in% names_high

### only need this if including lineage as covariate ###
#sample.info <- load.from.taiga(data.name='public-19q4-93d9', data.version=21, data.file='sample_info')
#rownames(sample.info) <- sample.info$DepMap_ID #if using CCLE names

source("running_limma.R")
lim_res <- run_lm_stats_limma(sub, in_group)#,covar = sample.info[rownames(sub),'lineage',drop=F])

# save output
write.table(lim_res,file='/Users/mkazachk/Documents/Uris_paper/Limma/drug_comparisons/prism_NEW.csv',sep = ",",row.names = FALSE)
