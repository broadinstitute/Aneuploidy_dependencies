if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

options(repos = c(
  "https://iwww.broadinstitute.org/~datasci/R-packages",
  "https://cran.cnr.berkeley.edu"))
install.packages('cdsr')

library(taigr)

ctd.matrix <- load.from.taiga(data.name='aneuploidy-data-d0b9', data.version=7, data.file='ctd_matrix')
aneuploidy.data <- load.from.taiga(data.name='aneuploidy-data-d0b9', data.version=18, data.file='aneuploidy_data_NEW_COMPACT')

high_aneuploidy = aneuploidy.data[,'many_arm_events'] == TRUE
low_aneuploidy = aneuploidy.data[,'many_arm_events'] == FALSE

# get intersecting lines
both = intersect(rownames(ctd.matrix),aneuploidy.data[,'CCLE_ID'])
sub = ctd.matrix [both,]

names_high = aneuploidy.data[high_aneuploidy,'CCLE_ID']
in_group = row.names(sub) %in% names_high

### only need this if including lineage as covariate ###
#sample.info <- load.from.taiga(data.name='public-19q4-93d9', data.version=21, data.file='sample_info')
#rownames(sample.info) <- sample.info$CCLE_Name #if using CCLE names

source("running_limma.R")
lim_res <- run_lm_stats_limma(sub, in_group)#,covar = sample.info[rownames(sub),'lineage',drop=F])

# write results to file
write.table(lim_res,file='/Users/mkazachk/Documents/Uris_paper/Limma/drug_comparisons/CTD_NEW.csv',sep = ",",row.names = FALSE)
