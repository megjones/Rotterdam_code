library(minfi)
library(sva)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library("lumi")
library("methylumi")
library(RCurl)
setwd("~/ewas3rdround")
source("deconvolution-utils.R") ## sometimes this gives an error - if so just open it and run the whole thing


load("sorted_betas_BMIQ_combat_together.rdata")
 load("sorted_pdat_BMIQ_combat_together.rdata")
 load("WB_pdat_BMIQ_combat_together.rdata")
 load("WB_betas_BMIQ_comabt_alloutliersremoved.rdata")
 
 ## remove samples with missing FACS
 validation_betas.combat<-  validation_betas.combat[,!colnames(validation_betas.combat)%in%c("7668610116_R03C01", "7668610134_R02C01", "7668610148_R02C02")]

 
 beta<-as.data.frame(validation_betas.combat)
 
 x <- getURL("https://raw.githubusercontent.com/redgar598/Tissue_Invariable_450K_CpGs/master/Invariant_Blood_CpGs.csv")
 y <- read.csv(text = x)
 beta_invariable<-beta[which(rownames(beta)%in%y$CpG),]#110228/114204 of the independent invariable sites are in beta
 Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
 beta_ref_range<-sapply(1:nrow(beta_invariable), function(x) Variation(beta_invariable[x,]))
 Invariable_in_beta<-beta_invariable[which(beta_ref_range<0.05),]
 invar_in_beta_and_independent<-intersect(y$CpG, rownames(Invariable_in_beta)) #106613/114204 (96.1%)
 
 save(invar_in_beta_and_independent, file="invariable_cordblood_CpGs.Rdata")

 
 validation_betas.combat<- validation_betas.combat[which(!(rownames(validation_betas.combat) %in% invar_in_beta_and_independent)),]
 save(validation_betas.combat, file="WB_betas_BMIQ_comabt_alloutliersremoved.rdata")
msorted_betas.combat<- msorted_betas.combat[which(!(rownames(msorted_betas.combat) %in% invar_in_beta_and_independent)),]


pd.sort<- msorted_pd.combat

CELL_TYPES = c("Gran", "Mono", "Bcell", "CD4T", "CD8T", "nRBC", "NK")


pd.sort$CellType <- to_canonical_cell_types(v = pd.sort$CellType,
                                            cell_types = c("WB", "CBMC", "CD4T", "CD8T", "G", "Mo", "B", "NK", "RBC"),
                                            canon_cell_types = c("WB", "CBMC", "CD4T", "CD8T", "Gran", "Mono", "Bcell", "NK", "nRBC"))

dm_probes_by_ct <- select_probes_t_separation_by_ct(betas=msorted_betas.combat,
                                                    cell_type_ind=pd.sort$CellType,
                                                    cell_types=CELL_TYPES, N=100)

probes <- unlist(dm_probes_by_ct)
save(probes, file="cord_cell_type_specific_probes.rdata")
library(tidyr)
cord_profile <- construct_profiles(betas=msorted_betas.combat,
                                   probes=probes,
                                   cell_types=CELL_TYPES, 
                                   cell_type_ind=pd.sort$CellType)
save(cord_profile, file="cord_profile.rdata")

est_cell_counts <- minfi:::projectCellType(validation_betas.combat[rownames(cord_profile),],
                                           cord_profile,
                                           nonnegative = T,
                                           lessThanOne = T)


save(est_cell_counts, file="GenR_cord_deconvolution_predicted_celltypes.rdata")


### Output tstats for all cell types and all probes
probe_coefs<- unlist(dm_probes_by_ct)

output_tstats_by_celltype <- function(betas, target_ct, cell_type_ind,
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
  output<- list(tstats)
}

# Wrapper for select_probes_t
output_tstats <- function(betas, cell_type_ind, cell_types) {
  ret <- lapply(cell_types, function(ct){output_tstats_by_celltype(betas=betas, 
                                                                   target_ct=ct,
                                                                   cell_type_ind=cell_type_ind,
                                                                   cell_types=cell_types )})
  
  names(ret) <- cell_types
  ret
}

probes_tstats <- output_tstats(betas=msorted_betas.combat,
                               cell_type_ind=pd.sort$CellType,
                               cell_types=CELL_TYPES)

save(probes_tstats, file="Cord_blood_celltype_ttest_statistics.Rdata")



