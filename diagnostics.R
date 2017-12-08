run_diagnostics <- function(betas1, pd1, betas2, pd2, betas12, pd12){
  # Distance Summary
  run_difference(betas1, pd1, betas2, pd2, betas12, pd12)
}

## Diff Summary
run_difference <- function(betas1, pd1, betas2, pd2){
  probes <- intersect(rownames(betas1), rownames(betas2))
  diff_summary <- diag_difference_summary(betas1[probes,], pd1$CellType,
                                          betas2[probes,], pd2$CellType,
                                          include_probes=F)
  
  # Filter out the null bucket.
  nbuckets <- length(levels(diff_summary$Interval))
  mid_idx <- (nbuckets %/% 2) + 1
  no_change_level <- levels(diff_summary$Interval)[mid_idx]
  
  plot <- diff_summary %>% filter(Interval != no_change_level) %>%
    ggplot(aes(Interval, Prop)) +
    geom_bar(stat="identity") +
    facet_grid( ~ CellType) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
    print(plot)
}

## Difference in mean betas per probe, summarized into buckets.
diag_difference_summary <- function(betas1, cell_types1, betas2, cell_types2,
                                    nbuckets=9, include_probes=F) {
  shared_cts <- intersect(cell_types1, cell_types2)
  res <- data.frame()
  probe_buckets <- list()
  for(ct in shared_cts){
    idx1 <- cell_types1 == ct
    b1 <- betas1[,idx1]
    idx2 <- cell_types2 == ct
    b2 <- betas2[,idx2]
    
    mean_diffs <- .probe_mean_difference(b1, b2, nbuckets)
    
    if(include_probes){
      probe_buckets[[ct]] <- mean_diffs
    }
    temp <- data.frame(table(mean_diffs))
    names(temp) <- c("Interval", "Count")
    temp$CellType = ct
    res <- rbind(res, temp)
  }
  
  res <- res %>% group_by(CellType) %>% mutate(Prop=Count/sum(Count))
  
  if(include_probes){
    probe_buckets
  } else{
    res
  }
}

.probe_mean_difference <- function(betas1, betas2, nbuckets) {
  # Check that rownames match.
  stopifnot(rownames(betas1) == rownames(betas2))
  
  probe_means1 <- rowMeans(betas1)
  probe_means2 <- rowMeans(betas2)
  diffs <- probe_means1 - probe_means2
  bins <- cut(diffs, seq(from=-.5, to=.5, length.out=(nbuckets + 1)))
  # Put the probe names back.
  names(bins) <- names(diffs)
  # Remove NA and NaN 
  bins <- bins[is.finite(bins)]
  
  bins
}

## PCA - From Kobor Lab
run_pca <- function(betas12, pd12, cat_vars=c("Slide", "Sex", "Row", "Col")){
  for(ct in ADULT_CELL_TYPES){
    ct_idx <- pd12$CellType == ct 
    complete <- complete.cases(betas12)
    PCA_full<-princomp(betas12[complete,ct_idx]) # only using the 417619 complete rows for PCA
    Loadings<-as.data.frame(unclass(PCA_full$loadings))
    vars <- PCA_full$sdev^2
    Importance<-vars/sum(vars)
    heat_scree_plot(Loadings, Importance,
                    title=ct,
                    meta_categorical = pd12[ct_idx, cat_vars],
                    meta_continuous = data.frame(RNum=runif(n = nrow(pd12[ct_idx,])))) %>%
      print()
  }
}
  
heat_scree_plot <- function(Loadings, Importance,
                            PCs_to_view=5, meta_categorical, meta_continuous,
                            title){
  ord <- c(seq(1:sum(ncol(meta_categorical), ncol(meta_continuous))))
  
  adjust<-1-Importance[1]
  pca_adjusted<-Importance[2:length(Importance)]/adjust
  pca_df<-data.frame(adjusted_variance=pca_adjusted, PC=seq(1:length(pca_adjusted)))
  
  scree<-ggplot(pca_df[which(pca_df$PC<(PCs_to_view+1)),],aes(PC,adjusted_variance))+geom_bar(stat = "identity",color="black",fill="grey")+theme_bw()+
    theme(axis.text = element_text(size =12),
          axis.title = element_text(size =15),
          plot.margin=unit(c(1.25,1.6,0.2,3),"cm"))+ylab("Adjusted Variance")+
    scale_x_continuous(breaks = seq(1,PCs_to_view,1)) + ggtitle(title)
  
  
  #### Heat
  ## correlate meta with PCS
  ## Run anova of each PC on each meta data variable
  
  
  aov_PC_meta <- lapply(1:ncol(meta_categorical), function(covar) sapply(1:ncol(Loadings), 
                                                                         function(PC) summary(aov(Loadings[, PC] ~ meta_categorical[, covar]))[[1]]$"Pr(>F)"[1]))
  cor_PC_meta <- lapply(1:ncol(meta_continuous), function(covar) sapply(1:ncol(Loadings), 
                                                                        function(PC) (cor.test(Loadings[, PC], as.numeric(meta_continuous[, 
                                                                                                                                          covar]), alternative = "two.sided", method = "spearman", na.action = na.omit)$p.value)))
  names(aov_PC_meta) <- colnames(meta_categorical)
  names(cor_PC_meta) <- colnames(meta_continuous)
  aov_PC_meta <- do.call(rbind, aov_PC_meta)
  cor_PC_meta <- do.call(rbind, cor_PC_meta)
  aov_PC_meta <- rbind(aov_PC_meta, cor_PC_meta)
  aov_PC_meta <- as.data.frame(aov_PC_meta)
  
  #adjust
  aov_PC_meta_adjust<-aov_PC_meta[,2:ncol(aov_PC_meta)]
  
  #reshape
  avo<-aov_PC_meta_adjust[,1:PCs_to_view]
  avo_heat_num<-apply(avo,2, as.numeric)
  avo_heat<-as.data.frame(avo_heat_num)
  avo_heat$meta<-rownames(avo)
  avo_heat_melt<-melt(avo_heat, id=c("meta"))
  
  # cluster meta data
  meta_var_order<-unique(avo_heat_melt$meta)[rev(ord)]
  avo_heat_melt$meta <- factor(avo_heat_melt$meta, levels = meta_var_order)
  
  # color if sig
  avo_heat_melt$Pvalue<-sapply(1:nrow(avo_heat_melt), function(x) if(avo_heat_melt$value[x]<=0.001){"<=0.001"}else{
    if(avo_heat_melt$value[x]<=0.01){"<=0.01"}else{
      if(avo_heat_melt$value[x]<=0.05){"<=0.05"}else{">0.05"}}})
  
  levels(avo_heat_melt$variable)<-sapply(1:PCs_to_view, function(x) paste("PC",x, sep="" ))
  
  heat<-ggplot(avo_heat_melt, aes(variable,meta, fill = Pvalue)) +
    geom_tile(color = "black",size=0.5) +
    theme_gray(8)+scale_fill_manual(values=c("#084594","#4292c6","#9ecae1","#deebf7"))+
    theme(axis.text = element_text(size =10, color="black"),
          axis.text.x = element_text(),
          axis.title = element_text(size =15),
          legend.text = element_text(size =14),
          legend.title = element_text(size =12),
          legend.position = c(1, 0.4), legend.justification = c(1,1),
          plot.margin=unit(c(0,2.25,1,1),"cm"))+
    xlab("Adjusted Principle Component")+ylab(NULL)
  
  grid.arrange(scree, heat, ncol=1, widths = 4, heights = c(2, 4))
}

## Hclust
run_hclust <- function(betas12, pd12){
  ## Hierarchical clustering.
  D <- dist(t(betas12))
  D.numeric <- as.matrix(D)
  
  col_colours <- cbind(CellType=as.numeric(factor(pd12$CellType)),
                       Sex=as.numeric(factor(pd12$Sex)),
                       Study=as.numeric(factor(pd12$Study)))
  print(heatmap3(D.numeric,
                 ColSideColors = col_colours,
                 main="Study"))
  
  col_colours <- cbind(CellType=as.numeric(factor(pd12$CellType)),
                       ChipRow=as.numeric(factor(pd12$Row)),
                       ChipCol=as.numeric(factor(pd12$Col)),
                       Chip=as.numeric(factor(pd12$Slide)))
  print(heatmap3(D.numeric, ColSideColors = col_colours, main="Chip"))
}

# MDS
run_mds <- function(betas12, pd12) {
  mdsPlot(betas12, sampGroups=pd12$CellType, main="CellType",
          pal=brewer.pal(11, "Spectral"), legendNCol = 2, legendPos="topright")
  mdsPlot(betas12,sampGroups=pd12$Study, main="Study",
          legendNCol = 2, legendPos="topright")
}

# T-test
run_t <- function(betas12, pd12){
  par(mfrow=c(2,3))
  for(ct in unique(pd12$CellType)){
    # Run a t-test coding the two studies against each other.
    ct_idx <- pd12$CellType == ct
    ct_betas <- betas12[,ct_idx]
    ct_pd <- pd12[ct_idx,]
    ct_m <- logit2(ct_betas)
    res <- rowttests(ct_m, fac=factor(ct_pd$Study))
    res$p.value %>% na.omit() %>% density() %>% plot(main=ct)
    summary(res$p.value)
  }
  par(mfrow=c(1,1))
}
