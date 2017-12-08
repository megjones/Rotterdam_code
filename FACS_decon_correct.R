setwd("~/ewas3rdround")

load("WB_betas_BMIQ_comabt_alloutliersremoved.rdata")
load("GenR_cord_deconvolution_predicted_celltypes.rdata")
load("genr_facs_counts.rdata")

facs<- validation_counts
facs_pcs<- prcomp(validation_counts, center=T, scale=T)$x
decon<- est_cell_counts
decon_pcs<- prcomp(decon, center=T, scale=T)$x

save(facs_pcs, file="facs_pcs.rdata")
save(decon_pcs, file="decon_pcs.rdata")
betas<- validation_betas.combat

for(a in c("facs", "facs_pcs", "decon", "decon_pcs")){
  diff<- as.data.frame(get(a))
  avebeta.lm<-apply(betas, 1, function(x){
    diff[colnames(betas),]->blood
    lm(as.formula(paste("x", "~",
                        paste(colnames(diff)[1:5], collapse = "+"),
                        sep = ""
    )),data=blood)
  })
  residuals<-t(sapply(avebeta.lm, function(x)residuals(summary(x))))
  colnames(residuals)<-colnames(betas)
  rownames(residuals)<- rownames(betas)
  adj.residuals<-residuals+matrix(apply(betas, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals)) ## figure out how to change the object name each time
  
  save(adj.residuals, file=paste(a, "_corrected_betas.Rdata", sep=""))}
