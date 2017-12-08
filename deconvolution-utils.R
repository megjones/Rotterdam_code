library(ggplot2)
library(dplyr)
library(tidyr)
library(genefilter)
library(quadprog)
library(BiocGenerics)

# TODO - If the type doesn't show up, just leave it alone!
to_canonical_cell_types <- function(v, cell_types, canon_cell_types=c("CD8T","CD4T","NK","Bcell","Mono", "Gran", "nRBC")){
  # Takes a vector of cell_types and an aligned canon_cell_types vector.
  # Translates v, a vector of cell types to canon_cell_types.
  idx <- match(v, cell_types)
  canon_cell_types[idx]
}

# Returns the plot for deconvolution accuracy with our generic datatype
# deconvolution_results is of type DataFrame
#   Fields:
#     - SampleName, CellType, True.Prop, Est.Prop
plot_deconv_accuracy <- function(deconvolution_results,
                                 plot_scales="fixed",
                                 summaries=c("pearson", "spearman", "mad")) {
  p <- ggplot(deconvolution_results, aes(True.Prop, Est.Prop)) +
       geom_point(aes(colour=SampleName)) +
       geom_abline(slope=1, intercept=0)
  
  # Generate the error measures and create labels out of them.
  errors <- deconvolution_results %>% 
              group_by(CellType) %>%
              summarize(mad = round(mean(abs(Est.Prop - True.Prop)), 2),
                        pearson=round(cor(Est.Prop, True.Prop, method="pearson"), 2),
                        spearman=round(cor(Est.Prop, True.Prop, method="spearman"), 2))
  
  # Create the labeller. Have statistic as a translation of long to short name.
  statistic = c(pearson="R", spearman="Rho", mad="MAD")
  labels <- errors$CellType
  for(s in summaries){
    summary <- paste(statistic[s], "=", errors[[s]], sep="")
    labels <- paste(labels, summary, sep="\n")
  }
  
  names(labels) <- errors$CellType
  p <- p + facet_wrap( ~ CellType, labeller=labeller(CellType = labels), scales=plot_scales) 
  p
}

# Helper function used in minfi
splitit <- function(x) {
  split(seq(along = x), x)
}

# Select probes based on their absolute separation.
select_probes_t_separation <- function(betas, target_ct, cell_type_ind,
                                       sort_abs=FALSE, N=50, p.value=1e-08, min.delta.beta=0,
                                       cell_types=c("Bcell", "Mono", "Gran", "CD4T", "CD8T", "NK", "nRBC")) {
  # Remove cell types we arn't interested in.
  keep <- cell_type_ind %in% cell_types
  betas <- betas[,keep]
  cell_type_ind <- cell_type_ind[keep]
  
  # Remove NA's
  betas <- na.omit(betas)
  
  tIndexes <- splitit(cell_type_ind)
  # Indicator var for two-group T-test
  target_idx <- tIndexes[[target_ct]]
  target_ind <- rep(0, ncol(betas))
  target_ind[target_idx] <- 1
  tstats <- genefilter::rowttests(betas, factor(target_ind))
  y <- tstats[tstats[, "p.value"] < p.value, ]
  y <- y[abs(y[,"dm"]) > min.delta.beta, ]
  y <- y[order(abs(y[,"dm"]), decreasing = T),]
  
  n <- min(nrow(y), N)
  ret <- row.names(y[1:n,])
  
  ret
}

# Wrapper for select_probes_t_separation.
select_probes_t_separation_by_ct <- function(betas, cell_type_ind, cell_types, N) {
  ret <- lapply(cell_types, function(ct){select_probes_t_separation(betas=betas, 
                                                                    target_ct=ct,
                                                                    cell_type_ind=cell_type_ind,
                                                                    cell_types=cell_types,
                                                                    N=N)})
  
  names(ret) <- cell_types
  ret
}

# Takes the betas, and gets the mean profile for each cell_type 
# in cell_types. Columns of betas should be labeled by cell_type_ind
.construct_profile <- function(betas, probes, cell_type, cell_type_ind) {
  ct_indexes <- splitit(cell_type_ind)
  ct_idx <- ct_indexes[[cell_type]]
  # Explicitly don't drop the dimension because it's possible to have only 1 sample of that ct
  ct_samples <- betas[probes,ct_idx,drop=FALSE]
  ct_ref <- as.matrix(rowMeans(ct_samples))
  colnames(ct_ref) <- cell_type
  
  ct_ref
}

construct_profiles <- function(betas, probes, cell_types, cell_type_ind) {
  res <- sapply(cell_types, function(ct){.construct_profile(betas=betas,
                                                    probes=probes,
                                                    cell_type=ct,
                                                    cell_type_ind=cell_type_ind)})
  
  row.names(res) <- probes
  res
}

# estimated_cell_counts is output for projectCellType
# validation_counts is tidy.
compare_estimation_to_truth <- function(estimated_cell_counts, validation_counts) {
  estimated_cell_counts <- as.data.frame(estimated_cell_counts)
  estimated_cell_counts$SampleName <- row.names(estimated_cell_counts)
  row.names(estimated_cell_counts) <- NULL
  
  data <- gather(estimated_cell_counts, CellType, Est.Prop, -SampleName) %>%
    arrange(SampleName)
  inner_join(data,
             validation_counts,
             by=c("SampleName", "CellType")) %>% arrange(SampleName, CellType)
}