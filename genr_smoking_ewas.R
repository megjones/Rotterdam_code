setwd("~/ewas3rdround")

### first, scale smoking corrected data
load("adj.residuals_sva.sup.smoke.Rdata")
load("adj.residuals_sva.unsup.smoke.Rdata")
load("adj.residuals_sva.ruv.smoke.Rdata")

for(a in c( "adj.residuals.ruv.smoke", 
           "adj.residuals.sva.sup.smoke",
           "adj.residuals.sva.unsup.smoke"
)){
  newbetas<- get(a)
  newbetas[newbetas>1]<- max(newbetas[newbetas<1])
  newbetas[newbetas<0]<- min(newbetas[newbetas>0])
  save(newbetas, file=paste(a, "corrected_betas_scaled.Rdata", sep=""))
}

## load alll data
load("WB_betas_BMIQ_comabt_alloutliersremoved.rdata")
uncorrected<- validation_betas.combat
load("decon_pcscorrected_betas_scaled.Rdata")
decon_pcs<- newbetas
load("deconcorrected_betas_scaled.Rdata")
decon<- newbetas
load("facs_pcscorrected_betas_scaled.Rdata")
facs_pcs<- newbetas
load("facscorrected_betas_scaled.Rdata")
facs<- newbetas
load("adj.residuals.reffreecellmixcorrected_betas_scaled.Rdata")
refactor<- newbetas
load("adj.residuals.refactorcorrected_betas_scaled.Rdata")
reffreecellmix <- newbetas
load("adj.residuals.ruv.smokecorrected_betas_scaled.Rdata")
ruv.smoke <- newbetas
load("adj.residuals.sva.sup.smokecorrected_betas_scaled.Rdata")
sva.sup.smoke <- newbetas
load("adj.residuals.sva.unsup.smokecorrected_betas_scaled.Rdata")
sva.unsup.smoke <- newbetas

metadata<- read.table("AnalysisfileUBC_including_smoking.dat", sep="\t", header=T)
rownames(metadata)<- metadata$Sample_ID
metadata<- metadata[!rownames(metadata)%in%c("8667045039_R02C01","7668610116_R03C01", "7668610134_R02C01", "7668610148_R02C02"),]


library(limma)
library(lumi)
mod.smoke<- model.matrix(~ SMOKE_SUSTAINED, data=metadata)

for(a in c("facs", 
           "facs_pcs",
           "decon", 
           "decon_pcs", 
           "refactor", 
           "reffreecellmix",
           "ruv.smoke", 
           "sva.sup.smoke", 
           "sva.unsup.smoke",
           "uncorrected"
)){
  bet<- get(a)
  mvals<- beta2m(bet)
  fit <- lmFit(mvals[,rownames(mod.smoke)], mod.smoke, na.action=na.omit)
  fit <- eBayes(fit)
  top <- topTable(fit, coef=c("SMOKE_SUSTAINED"), adjust="BH", number=Inf) # all significant hits
  db<-  rowMeans(bet[,metadata$SMOKE_SUSTAINED=="1"], na.rm=TRUE ) - rowMeans(bet[,metadata$SMOKE_SUSTAINED=="0"], na.rm=TRUE )
  top$db<- db
  top.filt.smoke<- top[,c("P.Value", "adj.P.Val", "db")]
  save(top.filt.smoke, file=paste(a, "_smoke_toptable.Rdata", sep=""))
}




## checking smoke toptables
load("facs_pcs_smoke_toptable.Rdata")
facs.pcs.top.smoke<- top.filt.smoke
load("facs_smoke_toptable.Rdata")
facs.top.smoke<- top.filt.smoke
load("decon_pcs_smoke_toptable.Rdata")
decon.pcs.top.smoke<- top.filt.smoke
load("decon_smoke_toptable.Rdata")
decon.top.smoke<- top.filt.smoke
load("refactor_smoke_toptable.Rdata")
refactor.top.smoke<- top.filt.smoke
load("reffreecellmix_smoke_toptable.Rdata")
reffreecellmix.top.smoke<- top.filt.smoke
load("ruv.smoke_smoke_toptable.Rdata")
ruv.top.smoke<- top.filt.smoke
load("sva.sup.smoke_smoke_toptable.Rdata")
sva.sup.top.smoke<- top.filt.smoke
load("sva.unsup.smoke_smoke_toptable.Rdata")
sva.unsup.top.smoke<- top.filt.smoke
load("uncorrected_smoke_toptable.Rdata")
uncorrected.top.smoke<- top.filt.smoke


####
smoke.top<- read.delim("smoke_toptable.txt", sep="\t", header=T, row.names = 1)
true_hits_smoke<-rownames(smoke.top)
length(true_hits_smoke)
toptables_list<- c("facs.pcs.top.smoke",
                   "facs.top.smoke", 
                   "decon.top.smoke", 
                   "decon.pcs.top.smoke", 
                   "refactor.top.smoke", 
                   "reffreecellmix.top.smoke", 
                   "ruv.top.smoke", 
                   "sva.sup.top.smoke", 
                   "sva.unsup.top.smoke",
                   "uncorrected.top.smoke"
)

### generate hits tables
hits.table.smoke<- data.frame(matrix(NA, 10, 6,
                                  dimnames=list(toptables_list,c("N_hits", 
                                                                 "True_Positives", 
                                                                 "False_Positives"))),
                           stringsAsFactors=F)


for(a in c("facs.pcs.top.smoke",
           "facs.top.smoke", 
           "decon.top.smoke", 
           "decon.pcs.top.smoke", 
           "refactor.top.smoke", 
           "reffreecellmix.top.smoke", 
           "ruv.top.smoke", 
           "sva.sup.top.smoke", 
           "sva.unsup.top.smoke",
           "uncorrected.top.smoke"
)){
  top<- get(a)
  hits.table.smoke[a, 1]<- length(intersect(rownames(top)[top$adj.P.Val< 0.05], 
                                         rownames(top)[abs(top$db)>0.02]))
  hits.table.smoke[a, 2]<- length(intersect(intersect(rownames(top)[top$adj.P.Val< 0.05], 
                                                   rownames(top)[abs(top$db)>0.02]),true_hits_smoke))
  hits.table.smoke[a,3]<- hits.table.smoke[a,1]-hits.table.smoke[a,2]
  }
hits.table.smoke$Sensitivity<- hits.table.smoke$True_Positives/(hits.table.smoke$True_Positives+hits.table.smoke$False_Negatives)*100
hits.table.smoke$Specificity<- (nrow(facs.pcs.top.smoke)-hits.table.smoke$False_Positives-hits.table.smoke$False_Negatives-hits.table.smoke$True_Positives)/
  ((nrow(facs.pcs.top.smoke)-hits.table.smoke$False_Positives-hits.table.smoke$False_Negatives-hits.table.smoke$True_Positives)+hits.table.smoke$False_Positives)*100


write.csv(hits.table.smoke, file="smoke_EWAS_hit_tables.csv")

