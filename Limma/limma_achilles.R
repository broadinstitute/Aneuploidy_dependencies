if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

options(repos = c(
  "https://iwww.broadinstitute.org/~datasci/R-packages",
  "https://cran.cnr.berkeley.edu"))
install.packages('cdsr')

library(taigr)

gene.effect <-  load.from.taiga(data.name='demeter2-achilles-5386', data.version=13, data.file='gene_effect')

#aneuploidy.data <- load.from.taiga(data.name='aneuploidy-data-d0b9', data.version=19, data.file='using_ten_percent_NEW')
#aneuploidy.data <- load.from.taiga(data.name='aneuploidy-data-d0b9', data.version=3, data.file='aneuploidy_data_compact')
#aneuploidy.data <- load.from.taiga(data.name='aneuploidy-data-d0b9', data.version=4, data.file='aneuploidy_data_compact_CCLE')
aneuploidy.data <- load.from.taiga(data.name='aneuploidy-data-d0b9', data.version=18, data.file='aneuploidy_data_NEW_COMPACT')

high_aneuploidy = aneuploidy.data[,'many_arm_events'] == TRUE
low_aneuploidy = aneuploidy.data[,'many_arm_events'] == FALSE

names_high = aneuploidy.data[high_aneuploidy,'CCLE_ID']
names_low = aneuploidy.data[low_aneuploidy,'CCLE_ID']

#get overlapping lines
both = intersect(rownames(gene.effect),aneuploidy.data[,'CCLE_ID'])

gene.effect_2 = gene.effect[both,]

# aneuploid group
in_group = row.names(gene.effect_2) %in% names_high


### only need this including lineage as covariate ###
#sample.info <- load.from.taiga(data.name='public-19q4-93d9', data.version=21, data.file='sample_info')
#rownames(sample.info) <- sample.info$CCLE_Name #if using CCLE names

source("running_limma.R")
lim_res <- run_lm_stats_limma(gene.effect_2, in_group)#,covar = sample.info[rownames(gene.effect_2),'lineage',drop=F])

# save the output
write.table(lim_res,file='lim_res_rnai_achilles_NEW.csv',sep = ",",row.names = FALSE)



