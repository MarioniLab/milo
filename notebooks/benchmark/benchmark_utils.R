### BENCHMARKING FUNCTIONS ###

library(SingleCellExperiment)
library(DAseq)
library(miloR)
library(tibble)
library(dplyr)
library(igraph)
# library(cydar)

## Set-up reticulate 4 MELD
reticulate::use_condaenv("emma_env", required=TRUE)
library(reticulate) ## development version of reticulate, or numba use breaks C stack

### SYNTHETIC LABELS ###

# Given a data_embedding, sample a simplex to weight each dimension to get a
# new PDF for each condition
# a: scaling coefficient in logit (lower a --> less extreme probabilities) 
.make_pdf <- function(data_embedding, a=0.2){
  # Create an array of values that sums to 1
  n_components = ncol(data_embedding)
  data_simplex = sort(runif(n = n_components-1))
  data_simplex = c(0, data_simplex, 1)
  data_simplex = diff(data_simplex)
  data_simplex = sample(data_simplex)
  # Weight each embedding component by the simplex weights
  sort_axis = rowSums(data_embedding * data_simplex)
  
  # Pass the weighted components through a logit
  pdf = 1/(1+exp(- a * sort_axis))
  if (sample(c(TRUE, FALSE), 1)){ pdf = 1 - pdf }
  return(pdf)  
}

## Smooth probabilities over KNN graph (to avoid having clusters with opposite sign DA)
.knn_smoothing <- function(graph, cond_probability, k=15, redDim='pca.corrected', d=10){
  # X_red_dim = reducedDim(sce, redDim)[,1:d]
  # graph = buildKNNGraph(t(X_red_dim), k = k)  
  adj = get.adjacency(graph)
  smooth_cond_probability <- (adj %*% cond_probability)/rowSums(adj)
  smooth_cond_probability
}

# Creates random differentially expressed regions over a dataset for benchmarking.
add_synthetic_labels <- function(sce, # SingleCellExperiment obj
                                 knn_graph, # for knn smoothing of probability values
                                 redDim='pca.corrected', # embedding to use to simulate differential abundance
                                 n_conditions=2, # number of conditions to simulate
                                 n_components=10, # number of components of embedding to use
                                 seed=42
){
  data_embedding = reducedDim(sce, redDim)[,1:n_components]
  set.seed(seed)
  # embedding data must be mean-centered
  data_embedding = t(scale(t(data_embedding), scale=FALSE))
  
  # Randomly flip sign of each embedding dimension
  data_embedding = apply(data_embedding, 2, function(x)  x*sample(c(-1, 1), 1) )
  
  conditions = paste0("Condition", 1:n_conditions)
  cond_probability = sapply(1:(length(conditions)-1), function(x) .make_pdf(data_embedding))
  
  # KNN Smoothing to avoid regions of the graph with opposite labels
  cond_probability = .knn_smoothing(knn_graph, cond_probability, redDim=redDim, d=n_components)
  
  # Normalize to sum to 1 for each cell
  # cond_probability <- t(apply(cond_probability, 1, function(x) x/sum(abs(x))))
  cond_probability = cbind(cond_probability, 1 - rowSums(cond_probability))
  colnames(cond_probability) = conditions
  
  # Generate labels for condition and replicates
  synth_labels <- sapply(1:nrow(cond_probability),  function(i) sample(colnames(cond_probability), size = 1, prob = cond_probability[i,]))
  synth_replicates <- rep(c("R1", "R2", "R3" ), 1000)
  synth_samples <- paste0(synth_labels, "_", synth_replicates)
  
  # Add synthetic labels and probabilities to colData
  colData(sce)[["synth_labels"]] <- synth_labels
  # colData(sce)[["synth_replicates"]] <- synth_replicates
  colData(sce)[["synth_samples"]] <- synth_samples
  colnames(cond_probability) <- paste0(colnames(cond_probability), "_prob")
  colData(sce)[colnames(cond_probability)] <- cond_probability
  return(sce)
}

### SYNTHETIC BATCH EFFECT ###

add_batch_effect <- function(embryo_sce, norm_sd=0.5){
  cellids_sample <- split(embryo_sce$cell, embryo_sce$synth_samples)
  X_pca <- reducedDim(embryo_sce, "pca.corrected")
  X_pca_batch <- X_pca

  for (b in names(cellids_sample)){
    batch_effect <- rnorm(ncol(X_pca), mean=0, sd = norm_sd)
    X_pca_batch[cellids_sample[[b]],] <- t(apply(X_pca_batch[cellids_sample[[b]],], 1, function(x) x + batch_effect))
  }
  
  reducedDim(embryo_sce, "pca_batch") <- X_pca_batch
  embryo_sce  
}

# ## Preprocessing/checking benchmarking input
# 
# prepare4bm <- function(){
#   ## Save library size for cydar
#   
#   ## Make sure there is a cell_id column in colData
#   
# }


### METHODS ###

## Milo

run_milo <- function(sce, condition_col, sample_col, reduced.dim="PCA",
                     k=15, d=30, prop=0.1, returnMilo = TRUE){
  ## Make design matrix
  design_df <- as.tibble(colData(sce)[c(sample_col, condition_col)]) %>%
    distinct() %>%
    column_to_rownames(sample_col)
  design <- formula(paste('~', condition_col, collapse = ' '))
  ## Build graph neighbourhoods
  milo <- Milo(sce)
  milo <- buildGraph(milo, k=k, d=d, reduced.dim = reduced.dim)
  milo <- makeNhoods(milo, prop = prop, k=k, d=d, reduced_dims = reduced.dim)
  ## Test DA
  milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample=sample_col)
  milo <- calcNhoodDistance(milo, d=d, reduced.dim = reduced.dim)
  DA_results <- testNhoods(milo, design = design, design.df = design_df)
  if (isTRUE(returnMilo)) {
    return(list(Milo=milo, DAres=DA_results))
  } else {
    DA_results 
  }
}

milo2output <- function(milo, da_res, out_type="continuous", alpha=0.1){
  if (out_type=="continuous") { 
    da.cell.mat <- milo@nhoods %*% da_res$logFC
    da.cell <- da.cell.mat[,1]
  } else {
    da.nhoods <- ifelse(da_res$SpatialFDR < alpha, ifelse(da_res$logFC > 0, "PosLFC", 'NegLFC'), "NotDA")
    da.nhoods.mat <- sapply(unique(da.nhoods), function(x) as.numeric(da.nhoods==x))
    da.cell.mat <- milo@nhoods %*% da.nhoods.mat
    da.cell <- apply(da.cell.mat, 1, function(x) colnames(da.cell.mat)[which.max(x)])
  }
  da.cell
}


## DAseq

run_daseq <- function(sce, k.vec, condition_col, cell_col="cell", reduced.dim = "PCA", d=30){
  condition_vec <- colData(sce)[[condition_col]]
  conds <- unique(condition_vec)
  cell.labels <- sce[[cell_col]]
  
  daseq_res <- getDAcells(X=reducedDim(sce, reduced.dim)[,1:d],
                          cell.labels=cell.labels,
                          labels.1=cell.labels[condition_vec == conds[2]],
                          labels.2=cell.labels[condition_vec == conds[1]],
                          k.vector=k.vec,
                          size=1)
  return(daseq_res)
  }

daseq2output <- function(sce, daseq_res, cell_col="cell", out_type="continuous"){
  if (out_type=="continuous") {
    da.cell <- daseq_res$da.pred
  } else {
    cell.labels <-  sce[[cell_col]]
    da.cell <- rep("NotDA", length(cell.labels))
    da.cell[daseq_res$da.up] <- "PosLFC"
    da.cell[daseq_res$da.down] <- "NegLFC"
  }
  da.cell
}

## MELD

run_meld_reticulate <- function(sce, condition_col, sample_col, reduced.dim="PCA", d=30, k=15) {
  ## Extract input from SingleCellExperiment
  X_red_dim = reducedDim(sce, reduced.dim)[,1:d]
  condition_vec <- colData(sce)[[condition_col]]
  sample_labels <- colData(sce)[[sample_col]]
  conditions <- sort(unique(condition_vec))
  ## Run MELD
  reticulate::source_python("./run_meld.py")
  py$run_meld(X_red_dim, sample_labels, conditions, k=k)
}

meld2output <- function(meld_res, likelihood_cutoff = 0.6, out_type="continuous"){
  conds <- colnames(meld_res)
  if (out_type=="continuous") {
    da.cell <- meld_res[conds[1]][,1]
  } else {
  da.cell <- rep("NotDA", nrow(meld_res))
  da.cell[meld_res[conds[1]] > likelihood_cutoff] <- "PosLFC"
  da.cell[meld_res[conds[2]] > likelihood_cutoff] <- "NegLFC"
  }
  da.cell
}

# ## Cydar
# 
# cd <- prepareCellData(processed.exprs)
# cd <- countCells(cd, tol=2, tol=2.0, filter=0, downsample=3)
# # do DA testing with edgeR
# cd.dge <- DGEList(assay(cd), lib.size=cd$totals)
# 
# # filter low abundance hyperspheres
# keep <- aveLogCPM(sim.dge) >= aveLogCPM(1, mean(sim.cydar$totals))
# sim.cydar <- sim.cydar[keep,]
# sim.dge <- sim.dge[keep,]
# 
# sim.design <- model.matrix(~Condition, data=test.meta[gsub(colnames(sim.cydar), pattern="\\.", replacement="_"), ])
# sim.dge <- estimateDisp(sim.dge, sim.design)
# sim.fit <- glmQLFit(sim.dge, sim.design)
# sim.res <- glmQLFTest(sim.fit, coef=2)
# 
# # control the spatial FDR
# cydar.res <- sim.res$table
# cydar.res$SpatialFDR <- spatialFDR(intensities(sim.cydar), sim.res$table$PValue)
# is.sig <- cydar.res$SpatialFDR <= 0.1
# summary(is.sig)

## Louvain clustering

run_louvain <- function(sce, condition_col, sample_col, k=15, d=30, reduced.dim="PCA"){
  ## Make design matrix
  design_df <- as.tibble(colData(sce)[c(sample_col, condition_col)]) %>%
    distinct() %>%
    column_to_rownames(sample_col)
  design <- formula(paste('~', condition_col, collapse = ' '))
  ## Louvain clustering
  X_red_dim = reducedDim(sce, reduced.dim)[,1:d]
  sce.graph <- buildKNNGraph(t(X_red_dim), k=k)
  louvain.clust <- cluster_louvain(sce.graph)
  louvain.clust.ids <- membership(louvain.clust)
  
  condition_vec <- colData(sce)[[condition_col]]
  sample_labels <- colData(sce)[[sample_col]]
  clust.df <- data.frame("cell_id"=colnames(sce), "Louvain.Clust"=as.character(louvain.clust.ids))
  clust.df$Sample <- sample_labels
  clust.df$Condition <- condition_vec
  
  louvain.model <- model.matrix(design, data=design_df)
  louvain.count <- as.matrix(table(clust.df$Louvain.Clust, clust.df$Sample))
  louvain.dge <- DGEList(counts=louvain.count, lib.size=log(colSums(louvain.count)))
  louvain.dge <- estimateDisp(louvain.dge, louvain.model)
  louvain.fit <- glmQLFit(louvain.dge, louvain.model, robust=TRUE)
  louvain.res <- as.data.frame(topTags(glmQLFTest(louvain.fit, coef=2), sort.by='none', n=Inf))
  
  clust.df$logFC <- louvain.res[clust.df$Louvain.Clust, 'logFC']
  clust.df$FDR <- louvain.res[clust.df$Louvain.Clust, 'FDR']
  return(clust.df)
  }

louvain2output <- function(louvain_res, out_type="continuous", alpha=0.1){
  if (out_type=="continuous") {
    da.cell <- louvain_res$logFC
  } else {
    da.cell <- ifelse(louvain_res$FDR < alpha, ifelse(louvain_res$logFC > 0, "PosLFC", 'NegLFC'), "NotDA")
  }
  da.cell
}

### RUN BENCHMARK ON SYNTHETIC LABELS ###
benchmark_da <- function(sce, condition_col='synth_labels', 
                         sample_col="synth_samples",
                         red_dim="pca.corrected",
                         params = list(milo = list(k=15),
                                       meld = list(k=15),
                                       daseq = list(k.vec=c(10,20,30, 40)),
                                       louvain = list(k=15)
                                       ),
                         d=30, out_type = "continuous"){
  ## Run milo
  milo_res <- run_milo(sce, condition_col=condition_col, sample_col=sample_col,reduced.dim = red_dim, d=d, k=params$milo$k)
  milo_out <- milo2output(milo_res$Milo, milo_res$DAres, out_type = out_type)
  ## Run DAseq
  daseq_res <- run_daseq(sce, k.vec=params$daseq$k.vec, condition_col, reduced.dim = red_dim, d=d)
  daseq_out <- daseq2output(sce, daseq_res, out_type = out_type)
  ## Run MELD
  meld_res <- run_meld_reticulate(sce, condition_col=condition_col, sample_col=sample_col,reduced.dim = red_dim, d=d, k=params$meld$k)
  meld_out <- meld2output(meld_res, out_type = out_type)
  ## Run louvain
  louvain_res <- run_louvain(sce, condition_col=condition_col, sample_col=sample_col,reduced.dim = red_dim, d=d, k=params$louvain$k)
  louvain_out <- louvain2output(louvain_res, out_type = out_type)
  ## Collect results + true labels
  bm <- data.frame(milo=milo_out, daseq=daseq_out, meld=meld_out, louvain=louvain_out)
  bm$true_prob <- sce$Condition2_prob 
  bm$true <- sce$true_labels
  long_bm <- pivot_longer(bm, cols = c(milo, daseq, meld, louvain), names_to='method', values_to="pred") 
  return(long_bm)
}

calculate_outcome <- function(long_bm){
  long_bm %>%
    mutate(outcome=case_when(true==pred & pred!="NotDA" ~ 'TP',
                             true!=pred & pred!="NotDA" ~ 'FP',
                             true!=pred & pred=="NotDA" ~ 'FN',
                             true==pred & pred=="NotDA"  ~ "TN"
    )) %>%
    group_by(method, outcome) %>%
    summarise(n=n()) %>%
    pivot_wider(id_cols=method, names_from=outcome, values_from=n, values_fill=0) %>%
    mutate(TPR=TP/(TP+FP), FPR=FP/(TP+FP), TNR=TN/(TN+FN), 
           Accuracy = (TP + TN)/(TP + TN + FP + FN),
           Recall = TP / (TP+FN)
    )
}


