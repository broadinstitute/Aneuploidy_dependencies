if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

options(repos = c(
  "https://iwww.broadinstitute.org/~datasci/R-packages",
  "https://cran.cnr.berkeley.edu"))
install.packages('cdsr')

library(taigr)

#for.lm.drive <- load.from.taiga(data.name='aneuploidy-data-d0b9', data.version=11, data.file='for_lm')
for.lm.achilles <- load.from.taiga(data.name='aneuploidy-data-d0b9', data.version=14, data.file='for_lm_achilles')
#gene.expression <- load.from.taiga(data.name='public-19q4-93d9', data.version=21, data.file='CCLE_expression')



# check if expression is a significant predictor of dependency when 
# accounting for lineage
summary_bub1b = lm(rnai_bub1b ~ BUB1B_exp + lineage , data = for.lm.achilles  )
summary_mad2l1 = lm(rnai_mad2l1 ~ MAD2L1_exp + lineage , data = for.lm.achilles  )


s_bub1b =  summary(summary_bub1b)[["coefficients"]]
s_mad2l1 =  summary(summary_mad2l1)[["coefficients"]]

write.table(s_bub1b,file='lm_bub1b_expression_achilles.csv',sep = ",",row.names = TRUE,col.names=TRUE)
write.table(s_mad2l1,file='lm_mad2l1_expression_achilles.csv',sep = ",",row.names = TRUE,col.names=TRUE)
