library(minfi)
library(sva)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library("lumi")
library("methylumi")
library(RCurl)

setwd("~/ewas3rdround")
source("deconvolution-utils_noimprovements.R") ## sometimes this gives an error - if so just open it and run the whole thing

load("sorted_pdat_BMIQ_combat_together.rdata")
load("genr_pdata.rdata")
load("genr_betas_noob_BMIQ.rdata")
load("2017.10.16_SCB_BMIQnorm.rdata")

pd.sort<- msorted_pd.combat
pd.test<-genr_pd.noob
msorted.betas<- betas(SCB.bm)
validation.betas<- genr.betas.noob.bmiq
over<- intersect(rownames(msorted.betas), rownames(validation.betas))
msorted.betas<- msorted.betas[over,]
validation.betas<- validation.betas[over,]


model.sort <- model.matrix(~ CellType, pd.sort)
msorted_betas.combat <- ComBat(msorted.betas, batch=pd.sort$Slide, mod = model.sort,par.prior = T, mean.only=FALSE)

validation_betas.combat <- ComBat(validation.betas, batch=pd.test$Array_number, mod = NULL,par.prior = T, mean.only=FALSE)


CELL_TYPES = c("Gran", "Mono", "Bcell", "CD4T", "CD8T", "nRBC", "NK")


pd.sort$CellType <- to_canonical_cell_types(v = pd.sort$CellType,
                                            cell_types = c("WB", "CBMC", "CD4T", "CD8T", "G", "Mo", "B", "NK", "RBC"),
                                            canon_cell_types = c("WB", "CBMC", "CD4T", "CD8T", "Gran", "Mono", "Bcell", "NK", "nRBC"))

dm_probes_by_ct <- select_probes_t_separation_by_ct(betas=msorted_betas.combat,
                                                    cell_type_ind=pd.sort$CellType,
                                                    cell_types=CELL_TYPES, N=50)

probes <- unlist(dm_probes_by_ct)
head(probes)
save(probes, file="cord_cell_type_specific_probes_noimprovements.rdata")
library(tidyr)
cord_profile <- construct_profiles(betas=msorted_betas.combat,
                                   probes=probes,
                                   cell_types=CELL_TYPES, 
                                   cell_type_ind=pd.sort$CellType)
save(cord_profile, file="cord_profile_noimprovements.rdata")


est_cell_counts <- minfi:::projectCellType(validation_betas.combat[rownames(cord_profile),],
                                           cord_profile,
                                           nonnegative = T,
                                           lessThanOne = T)


save(est_cell_counts, file="GenR_cord_deconvolution_predicted_celltypes_noimprovements.rdata")

