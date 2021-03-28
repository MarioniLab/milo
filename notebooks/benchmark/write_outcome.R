## Compute benchmark outcomes ##  
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(SingleCellExperiment)
  library(scran)
})

source('./benchmark_utils.R')

parser <- ArgumentParser()
parser$add_argument("benchmark_csv", type="character",
                    help = "path to csv file storing benchmark results")
parser$add_argument("p_thresh", type="double",
                    help = "probability thresh for true positives")
args <- parser$parse_args()

## Calculate outcome by setting a probability threshold for true DA (da_upper)
.outcome_by_prob <- function(benchmark_df, da_upper){
  DA_thresh <- 1 - da_upper
  benchmark_df <- benchmark_df %>%
    mutate(true = ifelse(true_prob < DA_thresh, "NegLFC", "NotDA")) 
  if (benchmark_df$method[1]=="meld") {
    benchmark_df <- mutate(benchmark_df, pred = ifelse(Condition2 < DA_thresh, "NegLFC", "NotDA")) 
  }
  ## Check if conditions were swapped in test
  pred_cor <- benchmark_df %>%
    mutate(pred=factor(pred, levels=c("NegLFC", "NotDA", "PosLFC"), ordered = TRUE)) %>%
    mutate(pred=as.numeric(pred)) %>%
    summarise(cor(pred, true_prob))
  ## Swap outcomes if so
  if (!is.na(pred_cor)) {
    if (pred_cor < -0.1) {
      benchmark_df <- mutate(benchmark_df, pred = ifelse(pred=="NegLFC", "PosLFC", ifelse(pred=="PosLFC", "NegLFC", "NotDA")))
    }
  }
  calculate_outcome(benchmark_df)
}

bm_file <- args$benchmark_csv
p_thresh <- args$p_thresh

## Make data frame w benchmark parameters
method <- str_split(bm_file, "\\.")[[1]][4]
enr <- str_remove_all(bm_file, ".+enr|_.+")
seed <- str_remove_all(bm_file, ".+seed|_.+")
batchEffect <- str_remove_all(bm_file, ".+batchEffect|_.+")
pop <- str_remove_all(bm_file, ".+pop_|_enr.+")
if (length(str_split(bm_file, "\\.")[[1]]) == 6) {
  k <- str_split(bm_file, "\\.")[[1]][5] %>% str_remove("k")
} else {
  k = 50
}

out_meta_df <- data.frame(file_id = str_remove_all(bm_file, ".+/benchmark_embryo_pop_|.csv")) %>%
  separate(col = file_id, sep = ".DAresults.", into=c("file_id", "method")) %>%
  separate(col = file_id, sep = "_enr", into=c("pop", "file_id")) %>%
  separate(col = file_id, sep = "_seed", into=c("enr", "file_id")) %>%
  mutate(seed=str_remove(seed, "seed"), batchEffect=str_remove_all(batchEffect, "batchEffect")) 

outcome_file <- bm_file %>% str_remove('.csv') %>% str_c(., ".outcome.p", p_thresh, ".csv")

## Save outcome 
read_csv(bm_file) %>%
  .outcome_by_prob(da_upper = p_thresh) %>%
  mutate(enr=enr, seed=seed, batchEffect=batchEffect, pop=pop, method=method, DA_thresh=p_thresh) %>%
  write.table(row.names = FALSE, sep = ",", col.names = FALSE, file=outcome_file)
