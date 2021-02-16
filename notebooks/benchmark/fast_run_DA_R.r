### Run DA methods in R ###

suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(SingleCellExperiment)
  library(scran)
})

source('./benchmark_utils.R')

parser <- ArgumentParser()
parser$add_argument("data_RDS", type="character",
                    help = "path to RDS storing SingleCellExperiment object")
parser$add_argument("coldata_obj", type="character",
                    help = "path to coldata")
parser$add_argument("pca_obj", type="character",
                    help = "path to pca")
parser$add_argument("method", type="character",
                    help = "DA method to use")
args <- parser$parse_args()

seed <- args$batch_seed
data_path <- args$data_RDS
coldata_path <- args$coldata_obj
pca_path <- args$pca_obj
DA_method <- args$method

## Load data
print("Loading dataset...")
sce <- readRDS(data_path)

## Load coldata and PCA
# outdir <- '/nfs/team205/ed6/data/milo_benchmark/synthetic_data/'
# outprefix <- str_c("benchmark_", data_id, "_pop_", pop, '_enr', pop_enr, "_seed", seed)
coldata <- read_csv(coldata_path) %>% column_to_rownames()
X_pca <-read_csv(pca_path) %>% column_to_rownames()

## Find DA probability x cell
k = 50
bm_params = list(
  milo = list(k=k),
  milo_batch = list(k=k),
  meld = list(k=k),
  daseq = list(k.vec=seq(k-(k/2), k+(k/2), by=10)),
  louvain = list(k=k)
  )

## Run DA method
out <- runDA(sce, X_pca, coldata = coldata, method = DA_method, params = bm_params)

## Save output 
bm_outdir <-'/nfs/team205/ed6/data/milo_benchmark/'
write_csv(out, str_c(bm_outdir, str_remove(pca_path, ".pca.csv"), ".DAresults.", DA_method, ".csv"))
