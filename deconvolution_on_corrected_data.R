## Deonconvolution on corrected data
library(minfi)
library(sva)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library("lumi")
library("methylumi")
library(ggridges)

setwd("~/ewas3rdround")


load("cord_profile.rdata")
load("decon_pcs_corrected_betas.Rdata")
decon_pcs<- adj.residuals
load("decon_corrected_betas.Rdata")
decon<- adj.residuals
load("facs_pcs_corrected_betas.Rdata")
facs_pcs<- adj.residuals
load("facs_corrected_betas.Rdata")
facs<- adj.residuals
load("adj.residuals_reffreecellmix.Rdata")
load("adj.residuals_refactor.Rdata")
load("adj.residuals_sva.sup.ga.Rdata")
load("adj.residuals_ruv.sex.Rdata")
load("adj.residuals_ruv.ga.Rdata")
load("adj.residuals_sva.unsup.ga.Rdata")
load("adj.residuals_sva.sup.sex.Rdata")
load("adj.residuals_sva.sup.sex.Rdata")
load("adj.residuals_sva.unsup.sex.Rdata")

for(a in c("facs", 
           "facs_pcs",
           "decon", 
           "decon_pcs", 
           "adj.residuals.refactor", 
           "adj.residuals.reffreecellmix",
           "adj.residuals.ruv.ga", 
           "adj.residuals.ruv.sex",
           "adj.residuals.sva.sup.ga",
           "adj.residuals.sva.sup.sex",
           "adj.residuals.sva.unsup.ga",
           "adj.residuals.sva.unsup.sex"
)){
  testbetas<- get(a)
  est_cell_counts <- minfi:::projectCellType(testbetas[rownames(cord_profile),],
                                             cord_profile,
                                             nonnegative = T,
                                             lessThanOne = T)
  save(est_cell_counts, file=paste(a, "_predicted_cell_types.Rdata", sep=""))
}



### test the resulting cell counts
load("adj.residuals.sva.unsup.sex_predicted_cell_types.Rdata")
sva.unsup.sex<- est_cell_counts
load("adj.residuals.sva.unsup.ga_predicted_cell_types.Rdata")
sva.unsup.ga<- est_cell_counts
load("adj.residuals.sva.sup.sex_predicted_cell_types.Rdata")
sva.sup.sex<- est_cell_counts
load("adj.residuals.sva.sup.ga_predicted_cell_types.Rdata")
sva.sup.ga<- est_cell_counts
load("adj.residuals.ruv.sex_predicted_cell_types.Rdata")
ruv.sex<- est_cell_counts
load("adj.residuals.ruv.ga_predicted_cell_types.Rdata")
ruv.ga<- est_cell_counts
load("adj.residuals.reffreecellmix_predicted_cell_types.Rdata")
cellmix<- est_cell_counts
load("adj.residuals.refactor_predicted_cell_types.Rdata")
refactor<- est_cell_counts
load("facs_predicted_cell_types.Rdata")
facs.counts<- est_cell_counts
load("facs_pcs_predicted_cell_types.Rdata")
facs.pcs.counts<- est_cell_counts
load("decon_predicted_cell_types.Rdata")
decon.counts<- est_cell_counts
load("decon_pcs_predicted_cell_types.Rdata")
decon.pcs.counts<- est_cell_counts

mads<- data.frame(FACS_counts=apply(facs.counts, 2, mad),
                  FACS_pcs=apply(facs.pcs.counts, 2, mad),
                  Decon_counts=apply(decon.counts, 2, mad),
                  Decon_pcs=apply(decon.pcs.counts, 2, mad), 
                  sva.unsup.sex= apply(sva.unsup.sex, 2, mad), 
                  sva.unsup.ga= apply(sva.unsup.ga, 2, mad), 
                  sva.sup.sex= apply(sva.sup.sex, 2, mad), 
                  sva.sup.ga= apply(sva.sup.ga, 2, mad), 
                  ruv.sex= apply(ruv.sex, 2, mad),
                  ruv.ga= apply(ruv.ga, 2, mad), 
                  cellmix= apply(cellmix, 2, mad), 
                  refactor=apply(refactor, 2, mad)
)

mads$celltypes<- rownames(mads)
save(mads, file="decon_oncorrected_mads.rdata")

