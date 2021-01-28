### Benchmarking methods with synthetic labels and batch effects on real data

suppressPackageStartupMessages({
    library(argparse)
    library(tidyverse)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
})

source('./benchmark_utils.R')

parser <- ArgumentParser()
parser$add_argument("data_RDS", type="character",
                    help = "path to RDS storing SingleCellExperiment object")
parser$add_argument("batch_seed", type="integer",
                    help = "Seed 4 batch effect")
parser$add_argument("--pop_enrichment", type="integer", default=0.7,
                    help = "Max condition probability in DA population")
parser$add_argument("--max_size", type="integer", default=1000,
                    help = "Min number of cells in population to select")
parser$add_argument("--k", type="integer", default=50,
                    help = "K parameter")
args <- parser$parse_args()

seed <- args$batch_seed
data_path <- args$data_RDS
k <- args$k
max_size <- args$max_size
pop_enr <- args$pop_enrichment

## Load data
print("Loading dataset...")
sce <- readRDS(data_path)

## Select population to simulate DA by size and save

sized_pops = names(table(sce$celltype))[table(sce$celltype) < max_size]
pop = sample(sized_pops, 1)

sce <- add_synthetic_labels_pop(sce, pop=pop, pop_column = "celltype", seed=seed, pop_enr=pop_enr)
true_labels <- ifelse(sce$Condition2_prob < 0.4, "NegLFC", ifelse(sce$Condition2_prob > 0.6, "PosLFC", "NotDA"))
colData(sce)[["true_labels"]] <- true_labels

## Simulate batch effects of different magnitude
set.seed(seed)
print("Simulating batch effects...")
bm_sce_ls <- lapply(c(0, 0.25, 0.75, 1), function(sd){
  sce_be <- add_batch_effect(sce, batch_col = "synth_batches", norm_sd=sd)
  sce_be$norm_sd <- sd
  sce_be
})

## Find DA probability x cell 
bm_params = list(
  milo = list(k=k),
  meld = list(k=k),
  daseq = list(k.vec=seq(k, k+50, by=10)),
  louvain = list(k=k)
  )

print("Saving params...")
outname_params <- str_c("/nfs/team205/ed6/data/milo_benchmark/benchmarkBatch_embryo_labels", seed, "_params.RDS")
saveRDS(bm_params, outname_params)

print("Running DA methods...")
i_bm = 1
for (x in bm_sce_ls){
  long_bm <- benchmark_da(x, out_type = "labels", red_dim = "pca_batch", d=30, params = bm_params)
  long_bm <- mutate(long_bm, batch_sd=x$norm_sd[1])
  long_bm <- mutate(long_bm, pop = pop, pop_enr=pop_enr, pop_size=sum(sce[["celltype"]]==pop) )
  outname <- str_c("/nfs/team205/ed6/data/milo_benchmark/benchmarkBatch_embryo_pop_size", max_size, '_enr', pop_enr, "_seed", seed, "_batchEffect", i_bm, ".csv")
  write_csv(long_bm, outname)
  i_bm = i_bm + 1 
}


