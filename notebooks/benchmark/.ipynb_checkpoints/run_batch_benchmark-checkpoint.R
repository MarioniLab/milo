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
parser$add_argument("labels_seed", type="integer",
                    help = "Seed 4 synthetic labels")
parser$add_argument("--k", type="integer", default=15,
                    help = "K parameter")
args <- parser$parse_args()

seed <- args$labels_seed
data_path <- args$data_RDS

## Load data
print("Loading dataset...")
sce <- readRDS(data_path)

## Add synthetic labels
print("Generating synthetic labels...")
# Build KNN graph for smoothing
X_red_dim = reducedDim(sce, "pca.corrected")[,1:30]
graph = buildKNNGraph(t(X_red_dim), k = 15)  

sce <- add_synthetic_labels(sce, n_components = 10,  redDim='pca.corrected', seed=seed, knn_graph = graph)
true_labels <- ifelse(sce$Condition2_prob < 0.4, "NegLFC", ifelse(sce$Condition2_prob > 0.6, "PosLFC", "NotDA"))
colData(sce)[["true_labels"]] <- true_labels

## Simulate batch effects of different magnitude
print("Simulating batch effects...")
bm_sce_ls <- lapply(c(0, 0.1, 0.3, 0.5, 0.7, 1), function(sd){
  sce_be <- add_batch_effect(sce, norm_sd=sd)
  sce_be$norm_sd <- sd
  sce_be
  })

## Find DA probability x cell 
bm_params = list(
  milo = list(k=15),
  meld = list(k=15),
  daseq = list(k.vec=c(10,20,30, 40)),
  louvain = list(k=15)
  )

print("Saving params...")
outname_params <- str_c("/nfs/team205/ed6/data/milo_benchmark/benchmarkBatch_embryo_labels", seed, "_params.RDS")
saveRDS(bm_params, outname_params)

print("Running DA methods...")
i_bm = 1
for (x in bm_sce_ls){
  long_bm <- benchmark_da(x, out_type = "labels", red_dim = "pca_batch", d=30, params = bm_params)
  long_bm <- mutate(long_bm, batch_sd=x$norm_sd[1])
  outname <- str_c("/nfs/team205/ed6/data/milo_benchmark/benchmarkBatch_embryo_labels", seed, "_batchEffect", i_bm, ".csv")
  write_csv(long_bm, outname)
  i_bm = i_bm + 1 
}


