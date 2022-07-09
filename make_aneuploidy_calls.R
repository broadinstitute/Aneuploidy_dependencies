library(magrittr)
library(tidyverse)
library(taigr)

out_dir <- '~/CPDS/data/CCLE/CCLE_arm_calls'

# Helper functions --------------------------------------------------------

#return boolean vector which is true if data is greater than upper threshold, and FALSE if lower than lower threshold
#use as_quantiles = T to specify whether the thresholds should be interpreted as quantiles
two_threshold_split <- function(data, thresh_vals = c(0.25, 0.75), as_quantiles = T) {
  group <- rep(NA, length(data))
  if (as_quantiles) {
    thresh_vals = quantile(data, thresh_vals, na.rm=T)
  }
  group[data >= thresh_vals[2]] <- TRUE
  group[data <= thresh_vals[1]] <- FALSE
  group
}

arm_calls <- read_csv(file.path(out_dir, 'CCLE_arm_call_matrix_withCCLE.csv')) %>% 
  dplyr::select(-CCLE_ID, -ploidy) %>% 
  tidyr::gather(key = 'chrom_arm', value = 'arm_call', -DepMap_ID)

# Compute aneuploidy features ---------------------------------------------
arm_CNVs_df <- arm_calls %>% 
  dplyr::mutate(arm_event = abs(arm_call)) %>% 
  dplyr::group_by(DepMap_ID) %>% 
  dplyr::summarise(num_arm_events = sum(arm_event, na.rm=T)) %>% 
  ungroup() %>% 
  dplyr::mutate(many_arm_events = two_threshold_split(num_arm_events),
                CCLE_ID = celllinemapr::arxspan.to.ccle(DepMap_ID))

#add in doubling time info
CCLE_annotations <- load.from.taiga(data.name='other-ccle2-c93e', data.version=2, data.file='Cell_lines_annotations_20181226')

RNAi_df <- read_csv('~/CPDS/data/Achilles/shRNA/TableS1_SampleInfo.csv') %>% 
  rename(RNAi_doubling_time = `Doubling time (hrs)`,
         CCLE_ID = Name) %>% 
  mutate(RNAi_doubling_time = str_replace_all(RNAi_doubling_time, '[> ]', ''),
         RNAi_doubling_time = str_replace_all(RNAi_doubling_time, 'hrs', ''),
         RNAi_doubling_time = as.numeric(RNAi_doubling_time))
CCLE_annotations %<>% left_join(RNAi_df, by = 'CCLE_ID')


arm_CNVs_df %<>% 
  left_join(CCLE_annotations %>% 
              dplyr::select(DepMap_ID, RNAi_doubling_time, CCLE_doubling_time = Doubling.Time.Calculated.hrs), 
            by = 'DepMap_ID')

write_csv(arm_CNVs_df, '~/CPDS/data/CCLE/aneuploidy_data.csv')


#merge in ploidy and WGD info and save out for taiga upload
CCLE_abs_table <- load.from.taiga(data.name='ccle-absolute-cn', data.version=5, data.file='CCLE_ABSOLUTE_combined_table')

arm_CNVs_df %>% 
    dplyr::select(DepMap_ID, `Aneuploidy score` = num_arm_events) %>% 
    left_join(CCLE_abs_table %>% 
                dplyr::select(DepMap_ID, Ploidy = ploidy, `Genome doublings`),
              by = 'DepMap_ID') %>% 
    write_csv('~/CPDS/data/CCLE/aneuploidy_data_for_taiga.csv')
