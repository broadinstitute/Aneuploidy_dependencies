library(magrittr)
library(tidyverse)
library(taigr)

out_dir <- '~/CPDS/data/CCLE/CCLE_arm_calls'

# Helper functions --------------------------------------------------------

#assigns chromosome arm to each segment using centromere range data
get_which_arm <- function(seg_start, seg_end, cent_start, cent_end) {
  output <- rep(NA, length(seg_start))
  seg_cent = 0.5*(seg_start + seg_end)
  output[seg_cent < cent_start] <- 'p'
  output[seg_cent > cent_end] <- 'q'
  return(output)
}

#call arm events based on comparing arm-level CN to ploidy values
arm_call <- function(CN_vals, ploidy_vals){
  arm_CN <- round(CN_vals)
  ploidy_CN <- round(ploidy_vals)
  arm_calls <- rep(0, times = length(CN_vals))
  arm_calls[arm_CN < ploidy_CN] <- -1
  arm_calls[arm_CN > ploidy_CN] <- 1
  arm_calls
}

#handle cases where a segment crosses over the centromere by splitting it there
split_cent_crosses <- function(seg_df) {
  cross_segs <- seg_df$Start < seg_df$centStart & seg_df$End > seg_df$centEnd
  left_side <- seg_df[cross_segs, ] %>% 
    mutate(End = 0.5*(centStart + centEnd))
  right_side <- seg_df[cross_segs, ] %>% 
    mutate(Start = 0.5*(centStart + centEnd))
  seg_df[cross_segs, ] <- left_side
  seg_df <- rbind(seg_df, right_side)
  return(seg_df)
}



# Load data ---------------------------------------------------------------

#load hg19 centromere info
cent_df <- rCGH::hg19 %>% 
  dplyr::mutate(chrom = as.character(chrom)) %>% 
  dplyr::select(chrom,
                centStart = centromerStart,
                centEnd = centromerEnd)

#load CCLE2 ABSOLUTE table and seg-file
CCLE_abs_table <- load.from.taiga(data.name='ccle-absolute-cn', data.version=5, data.file='CCLE_ABSOLUTE_combined_table')

seg_df <- load.from.taiga(data.name='ccle-absolute-cn', data.version=5, data.file='CCLE_ABSOLUTE_combined_segtab') %>% 
    dplyr::mutate(CN = Modal_HSCN_1 + Modal_HSCN_2) %>% 
  dplyr::select(chrom = Chromosome,
                Start,
                End,
                CN,
                CCLE_ID,
                DepMap_ID) %>% 
  dplyr::mutate(chrom = as.character(chrom)) %>% 
  left_join(cent_df, by = 'chrom') %>% 
  split_cent_crosses() %>% #split segments at centromeres
  dplyr::mutate(width = End-Start) #compute segment widths


# Calculate arm-calls -----------------------------------------------------

#assign arms to each segment and compute total arm lengths and relative weights
arm_calls <- seg_df %>% 
  dplyr::mutate(arm = get_which_arm(Start, End, centStart, centEnd)) %>% 
  dplyr::filter(!is.na(arm)) %>% 
  dplyr::group_by(DepMap_ID, chrom, arm) %>% 
  dplyr::summarize(wmed_CN = spatstat::weighted.median(CN, width, na.rm=T)) %>% 
  ungroup() %>% 
  left_join(
    CCLE_abs_table %>% dplyr::select(DepMap_ID, ploidy),
    by = 'DepMap_ID'
  ) %>% dplyr::mutate(arm_call = arm_call(wmed_CN, ploidy)) 


#save matrix of arm-calls
arm_ord <- c('1p', '1q', '2p', '2q', '3p', '3q', '4p', '4q', '5p', '5q', '6p', '6q', '7p', '7q', '8p', '8q', '9p', '9q', '10p', '10q', '11p', '11q', '12p', '12q', '13q', '14q', '15q', '16p', '16q', '17p', '17q', '18p', '18q', '19p', '19q', '20p', '20q', '21q', '22q')
arm_call_mat <- arm_calls %>% 
    dplyr::mutate(CCLE_ID = celllinemapr::arxspan.to.ccle(DepMap_ID),
                  chrom_arm = paste0(chrom, arm),
                  chrom_arm = factor(chrom_arm, levels = arm_ord)) %>% 
    dplyr::select(DepMap_ID, CCLE_ID, ploidy, chrom_arm, arm_call) %>% 
    tidyr::spread(key = 'chrom_arm', value = 'arm_call') 
write_csv(arm_call_mat, file.path(out_dir, 'CCLE_arm_call_matrix_withCCLE.csv'))
arm_call_mat %>% 
    dplyr::select(-CCLE_ID, -ploidy) %>% 
    write_csv(file.path(out_dir, 'CCLE_arm_call_matrix.csv'))

#write out files used from taiga 
write_csv(CCLE_abs_table, file.path(out_dir, 'CCLE_ABSOLUTE_combined_table.csv'))
write_csv(seg_df, file.path(out_dir, 'CCLE_ABSOLUTE_combined_segtab.csv'))
