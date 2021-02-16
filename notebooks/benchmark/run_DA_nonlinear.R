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
parser$add_argument("batch_seed", type="integer",
                    help = "Seed 4 batch effect")
parser$add_argument("population", type="character",
                    help = "Which cell type is DA?")
parser$add_argument("--pop_enrichment", type="double", default=0.7,
                    help = "Max condition probability in DA population")
parser$add_argument("--k", type="integer", default=50,
                    help = "K parameter")
args <- parser$parse_args()

seed <- args$batch_seed
data_path <- args$data_RDS
k <- args$k
pop <- args$population
pop_enr <- args$pop_enrichment

## Load data
print("Loading dataset...")
sce <- readRDS(data_path)

## Load coldata and PCA
outdir <- '/nfs/team205/ed6/data/milo_benchmark/synthetic_data/'
outprefix <- str_c("benchmark_embryo_pop_", pop, '_enr', pop_enr, "_seed", seed)
coldata <- read_csv(paste0(outdir, outprefix, ".coldata.csv")) %>% column_to_rownames()
X_pca <-read_csv(str_c(outdir, outprefix, "_batchEffectNonLinear30.pca.csv")) %>% column_to_rownames()

## Find DA probability x cell
bm_params = list(
  milo = list(k=k),
  milo_batch = list(k=k),
  meld = list(k=k),
  daseq = list(k.vec=seq(k-(k/2), k+(k/2), by=10)),
  louvain = list(k=k)
  )

## Run DA method
out <- runDA(sce, X_pca, coldata = coldata, method = "milo_batch", params = bm_params)

## Save output 
bm_outdir <-'/nfs/team205/ed6/data/milo_benchmark/'
write_csv(out, str_c(bm_outdir, outprefix, "_batchEffectNonLinear30.DAresults.milo_batch.csv"))
