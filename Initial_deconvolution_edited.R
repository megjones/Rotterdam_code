## Louie's deconvolution

#  Louie's deconvolution
library(minfi)
library(sva)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library("lumi")
library("methylumi")

source("~/UBC_files/deconvolution-utils.R") ## sometimes this gives an error - if so just open it and run the whole thing
source("~/UBC_files/diagnostics.R")


load("2017.10.16_SCB_BMIQnorm.rdata")
load("genr_pdata.rdata")
load("genr_betas_noob_BMIQ.rdata")

colnames(genr.betas.noob.bmiq) <-gsub("/home/040129/3roundEWAs-AllFiles/", "", colnames(genr.betas.noob.bmiq))
rownames(genr_pd.noob) <-gsub("/home/040129/3roundEWAs-AllFiles/", "", rownames(genr_pd.noob))
identical(rownames(genr_pd.noob), colnames(genr.betas.noob.bmiq))

msorted_betas<- betas(SCB.bm)
testbetas<- genr.betas.noob.bmiq
over<- intersect(rownames(msorted_betas), rownames(testbetas))
combined_betas <- cbind(msorted_betas[over,], testbetas[over,])

pd.sort<- pData(SCB.bm)
pd.test<-genr_pd.noob
combined_pd <- rbind(data.frame(Sample_Name=rownames(pd.sort), CellType=pd.sort$Tissue,  Array=pd.sort$Sentrix_Position, Slide=pd.sort$Sentrix_ID), 
          data.frame(Sample_Name=rownames(pd.test), CellType=rep("WB", nrow(pd.test)),  Array=pd.test$Position, Slide=pd.test$Array_number))
combined_pd$Study <- c(rep("Ref", times=ncol(msorted_betas)), rep("Val", ncol(testbetas)))

combined.m<- beta2m(combined_betas)
model <- model.matrix(~ CellType, combined_pd)
table(combined_pd$Slide, combined_pd$CellType)
res <- ComBat(combined.m, batch=combined_pd$Slide, mod = model,par.prior = T, mean.only=FALSE)

ref_idx <- combined_pd$Study == "Ref"
msorted_betas.combat <- m2beta(res[,ref_idx])## changed this from initial code to add m2beta part
msorted_pd.combat <- combined_pd[ref_idx,]
val_idx <- combined_pd$Study == "Val"
validation_betas.combat <- m2beta(res[,val_idx])## changed this from initial code to add m2beta part
validation_pd.combat <- combined_pd[val_idx,]
save(msorted_betas.combat, file="sorted_betas_BMIQ_combat_together.rdata")
save(msorted_pd.combat, file="sorted_pdat_BMIQ_combat_together.rdata")
save(validation_pd.combat, file="WB_pdat_BMIQ_combat_together.rdata")
save(validation_betas.combat, file="WB_betas_BMIQ_combat_together.rdata")

### repeated this bc have to remove invariable probes before deconvolution
load("sorted_betas_BMIQ_combat_together.rdata")
load("sorted_pdat_BMIQ_combat_together.rdata")
load("WB_pdat_BMIQ_combat_together.rdata")
load("WB_betas_BMIQ_combat_together.rdata")


### identify invariable probes in GenR
## identify invariable probes
beta<-as.data.frame(validation_betas.combat)

x <- getURL("https://raw.githubusercontent.com/redgar598/Tissue_Invariable_450K_CpGs/master/Invariant_Blood_CpGs.csv")
y <- read.csv(text = x)
beta_invariable<-beta[which(rownames(beta)%in%y$CpG),]#110228/114204 of the independent invariable sites are in beta
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
beta_ref_range<-sapply(1:nrow(beta_invariable), function(x) Variation(beta_invariable[x,]))
Invariable_in_beta<-beta_invariable[which(beta_ref_range<0.05),]
invar_in_beta_and_independent<-intersect(y$CpG, rownames(Invariable_in_beta)) #106613/114204 (96.1%)

save(invar_in_beta_and_independent, file="invariable_cordblood_CpGs.Rdata")


msorted_betas.combat<- msorted_betas.combat[!which(rownames(msorted_betas.combat %in% invar_in_beta_and_independent)),]
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



