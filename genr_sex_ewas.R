### First, rescale any values above 1 or below 0

load("/home/projects/genr/ewas3rdround/decon_pcs_corrected_betas.Rdata")
decon_pcs<- adj.residuals
load("/home/projects/genr/ewas3rdround/decon_corrected_betas.Rdata")
decon<- adj.residuals
load("/home/projects/genr/ewas3rdround/facs_pcs_corrected_betas.Rdata")
facs_pcs<- adj.residuals
load("/home/projects/genr/ewas3rdround/facs_corrected_betas.Rdata")
facs<- adj.residuals
load("/home/projects/genr/ewas3rdround/adj.residuals_reffreecellmix.Rdata")
load("/home/projects/genr/ewas3rdround/adj.residuals_refactor.Rdata")
load("/home/projects/genr/ewas3rdround/invariable_cordblood_CpGs.Rdata")
load("/home/projects/genr/ewas3rdround/adj.residuals_sva.sup.ga.Rdata")
load("/home/projects/genr/ewas3rdround/adj.residuals_ruv.sex.Rdata")
load("/home/projects/genr/ewas3rdround/adj.residuals_ruv.ga.Rdata")
load("/home/projects/genr/ewas3rdround/adj.residuals_sva.unsup.ga.Rdata")
load("/home/projects/genr/ewas3rdround/adj.residuals_sva.sup.sex.Rdata")
load("/home/projects/genr/ewas3rdround/adj.residuals_sva.sup.sex.Rdata")
load("/home/projects/genr/ewas3rdround/adj.residuals_sva.unsup.sex.Rdata")

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
  newbetas<- get(a)
  newbetas[newbetas>1]<- max(newbetas[newbetas<1])
  newbetas[newbetas<0]<- min(newbetas[newbetas>0])
  save(newbetas, file=paste(a, "corrected_betas_scaled.Rdata", sep=""))
}


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
load("adj.residuals.ruv.gacorrected_betas_scaled.Rdata")
ruv.ga <- newbetas
load("adj.residuals.ruv.sexcorrected_betas_scaled.Rdata")
ruv.sex <- newbetas
load("adj.residuals.sva.sup.gacorrected_betas_scaled.Rdata")
sva.sup.ga <- newbetas
load("adj.residuals.sva.sup.sexcorrected_betas_scaled.Rdata")
sva.sup.sex <- newbetas
load("adj.residuals.sva.unsup.gacorrected_betas_scaled.Rdata")
sva.unsup.ga <- newbetas
load("adj.residuals.sva.unsup.sexcorrected_betas_scaled.Rdata")
sva.unsup.sex <- newbetas


library(limma)
library(lumi)
load("genr_meta.rdata")
mod.sex<- model.matrix(~Sex, data=metadata)

for(a in c("facs", 
           "facs_pcs",
           "decon", 
           "decon_pcs", 
           "refactor", 
           "reffreecellmix",
           "ruv.sex", 
           "sva.sup.sex", 
           "sva.unsup.sex",
           "uncorrected"
)){
  bet<- get(a)
  mvals<- beta2m(bet)
  fit <- lmFit(mvals[,rownames(mod.sex)], mod.sex, na.action=na.omit)
  fit <- eBayes(fit)
  top <- topTable(fit, coef=c("Sex2"), adjust="BH", number=Inf) # all significant hits
  db<-  rowMeans(bet[,metadata$Sex=="2"], na.rm=TRUE ) - rowMeans(bet[,metadata$Sex=="1"], na.rm=TRUE )
  top$db<- db[rownames(top)]
  top.filt.sex<- top[,c("P.Value", "adj.P.Val", "db")]
  save(top.filt.sex, file=paste(a, "_sex_toptable.Rdata", sep=""))
}




## checking sex toptables
load("facs_pcs_sex_toptable.Rdata")
facs.pcs.top.sex<- top.filt.sex
load("facs_sex_toptable.Rdata")
facs.top.sex<- top.filt.sex
load("decon_pcs_sex_toptable.Rdata")
decon.pcs.top.sex<- top.filt.sex
load("decon_sex_toptable.Rdata")
decon.top.sex<- top.filt.sex
load("refactor_sex_toptable.Rdata")
refactor.top.sex<- top.filt.sex
load("reffreecellmix_sex_toptable.Rdata")
reffreecellmix.top.sex<- top.filt.sex
load("ruv.sex_sex_toptable.Rdata")
ruv.top.sex<- top.filt.sex
load("sva.sup.sex_sex_toptable.Rdata")
sva.sup.top.sex<- top.filt.sex
load("sva.unsup.sex_sex_toptable.Rdata")
sva.unsup.top.sex<- top.filt.sex
load("uncorrected_sex_toptable.Rdata")
uncorrected.top.sex<- top.filt.sex


####
load("sorted_sex_toptable.Rdata")
sig_hits_sex<- sort.ewas.sex.top2[intersect(rownames(sort.ewas.sex.top2)[sort.ewas.sex.top2$Pvalue< 0.000001], 
                                            rownames(sort.ewas.sex.top2)[abs(sort.ewas.sex.top2$db)>0.02]),]
true_hits_sex<-rownames(sig_hits_sex)
length(true_hits_sex)
toptables_list<- c("facs.pcs.top.sex",
                   "facs.top.sex", 
                   "decon.top.sex", 
                   "decon.pcs.top.sex", 
                   "refactor.top.sex", 
                   "reffreecellmix.top.sex", 
                   "ruv.top.sex", 
                   "sva.sup.top.sex", 
                   "sva.unsup.top.sex",
                   "uncorrected.top.sex"
)

### generate hits tables
hits.table.sex<- data.frame(matrix(NA, 10, 6,
                                   dimnames=list(toptables_list,c("N_hits", 
                                                                  "True_Positives", 
                                                                  "False_Positives", 
                                                                  "False_Negatives",
                                                                  "Spearman", 
                                                                  "Kendall"))),
                            stringsAsFactors=F)


for(a in c("facs.pcs.top.sex",
           "facs.top.sex", 
           "decon.top.sex", 
           "decon.pcs.top.sex", 
           "refactor.top.sex", 
           "reffreecellmix.top.sex", 
           "ruv.top.sex", 
           "sva.sup.top.sex", 
           "sva.unsup.top.sex",
           "uncorrected.top.sex"
)){
  top<- get(a)
  hits.table.sex[a, 1]<- length(intersect(rownames(top)[top$adj.P.Val< 0.05], 
                                          rownames(top)[abs(top$db)>0.02]))
  hits.table.sex[a, 2]<- length(intersect(intersect(rownames(top)[top$adj.P.Val< 0.05], 
                                                     rownames(top)[abs(top$db)>0.02]),true_hits_sex))
  hits.table.sex[a,3]<- hits.table.sex[a,1]-hits.table.sex[a,2]
  hits.table.sex[a,4]<- cor(sort.ewas.sex.top2$Pvalue, top[rownames(sort.ewas.sex.top2),"P.Value"], method="spearman", use="pairwise")
  hits.table.sex[a,5]<- cor(sort.ewas.sex.top2[1:1000,"Pvalue"], top[rownames(sort.ewas.sex.top2[1:1000,]),"P.Value"], method="kendall", use="pairwise")
}
hits.table.sex$Sensitivity<- hits.table.sex$True_Positives/(hits.table.sex$True_Positives+hits.table.sex$False_Negatives)*100
hits.table.sex$Specificity<- (nrow(facs.pcs.top.sex)-hits.table.sex$False_Positives-hits.table.sex$False_Negatives-hits.table.sex$True_Positives)/
  ((nrow(facs.pcs.top.sex)-hits.table.sex$False_Positives-hits.table.sex$False_Negatives-hits.table.sex$True_Positives)+hits.table.sex$False_Positives)*100


write.csv(hits.table.sex, file="Sex_EWAS_hit_tables.csv")

