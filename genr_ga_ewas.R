setwd("~/ewas3rdround")

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
load("adj.residuals.sva.sup.gacorrected_betas_scaled.Rdata")
sva.sup.ga <- newbetas
load("adj.residuals.sva.sup.gacorrected_betas_scaled.Rdata")
sva.sup.ga <- newbetas
load("adj.residuals.sva.unsup.gacorrected_betas_scaled.Rdata")
sva.unsup.ga <- newbetas
load("adj.residuals.sva.unsup.gacorrected_betas_scaled.Rdata")
sva.unsup.ga <- newbetas


library(limma)
library(lumi)
load("genr_meta.rdata")
mod.ga<- model.matrix(~GA, data=metadata)

for(a in c("facs", 
           "facs_pcs",
           "decon", 
           "decon_pcs", 
           "refactor", 
           "reffreecellmix",
           "ruv.ga", 
           "sva.sup.ga", 
           "sva.unsup.ga",
           "uncorrected"
)){
  bet<- get(a)
  mvals<- beta2m(bet)
  fit <- lmFit(mvals[,rownames(mod.ga)], mod.ga, na.action=na.omit)
  fit <- eBayes(fit)
  top <- topTable(fit, coef=c("GA"), adjust="BH", number=Inf) # all significant hits
  db<- sapply(1:nrow(bet), function(x) {
    z<-lm(unlist(bet[x,]) ~ metadata$GA)
    slope=z$coefficients[2]
    as.numeric(slope*(max(metadata$GA)-min(metadata$GA))) 
  })
  top$db<- db
  top.filt.ga<- top[,c("P.Value", "adj.P.Val", "db")]
  save(top.filt.ga, file=paste(a, "_ga_toptable.Rdata", sep=""))
}




## checking ga toptables
load("facs_pcs_ga_toptable.Rdata")
facs.pcs.top.ga<- top.filt.ga
load("facs_ga_toptable.Rdata")
facs.top.ga<- top.filt.ga
load("decon_pcs_ga_toptable.Rdata")
decon.pcs.top.ga<- top.filt.ga
load("decon_ga_toptable.Rdata")
decon.top.ga<- top.filt.ga
load("refactor_ga_toptable.Rdata")
refactor.top.ga<- top.filt.ga
load("reffreecellmix_ga_toptable.Rdata")
reffreecellmix.top.ga<- top.filt.ga
load("ruv.ga_ga_toptable.Rdata")
ruv.top.ga<- top.filt.ga
load("sva.sup.ga_ga_toptable.Rdata")
sva.sup.top.ga<- top.filt.ga
load("sva.unsup.ga_ga_toptable.Rdata")
sva.unsup.top.ga<- top.filt.ga
load("uncorrected_ga_toptable.Rdata")
uncorrected.top.ga<- top.filt.ga


####
load("sorted_ga_toptable.Rdata")
sig_hits_ga<- sort.ewas.ga.top[intersect(rownames(sort.ewas.ga.top)[sort.ewas.ga.top$Pvalue< 0.000001], 
                                            rownames(sort.ewas.ga.top)[abs(sort.ewas.ga.top$db)>0.02]),]
true_hits_ga<-rownames(sig_hits_ga)
length(true_hits_ga)
toptables_list<- c("facs.pcs.top.ga",
                   "facs.top.ga", 
                   "decon.top.ga", 
                   "decon.pcs.top.ga", 
                   "refactor.top.ga", 
                   "reffreecellmix.top.ga", 
                   "ruv.top.ga", 
                   "sva.sup.top.ga", 
                   "sva.unsup.top.ga",
                   "uncorrected.top.ga"
)

### generate hits tables
hits.table.ga<- data.frame(matrix(NA, 10, 6,
                                   dimnames=list(toptables_list,c("N_hits", 
                                                                  "True_Positives", 
                                                                  "False_Positives", 
                                                                  "False_Negatives",
                                                                  "Spearman", 
                                                                  "Kendall"))),
                            stringsAsFactors=F)


for(a in c("facs.pcs.top.ga",
           "facs.top.ga", 
           "decon.top.ga", 
           "decon.pcs.top.ga", 
           "refactor.top.ga", 
           "reffreecellmix.top.ga", 
           "ruv.top.ga", 
           "sva.sup.top.ga", 
           "sva.unsup.top.ga",
           "uncorrected.top.ga"
)){
  top<- get(a)
  hits.table.ga[a, 1]<- length(intersect(rownames(top)[top$adj.P.Val< 0.05], 
                                          rownames(top)[abs(top$db)>0.02]))
  hits.table.ga[a, 2]<- length(intersect(intersect(rownames(top)[top$adj.P.Val< 0.05], 
                                                    rownames(top)[abs(top$db)>0.02]),true_hits_ga))
  hits.table.ga[a,3]<- hits.table.ga[a,1]-hits.table.ga[a,2]
  hits.table.ga[a,4]<- cor(sort.ewas.ga.top2$Pvalue, top[rownames(sort.ewas.ga.top2),"P.Value"], method="spearman", use="pairwise")
  hits.table.ga[a,5]<- cor(sort.ewas.ga.top2[1:1000,"Pvalue"], top[rownames(sort.ewas.ga.top2[1:1000,]),"P.Value"], method="kendall", use="pairwise")
}
hits.table.ga$Sensitivity<- hits.table.ga$True_Positives/(hits.table.ga$True_Positives+hits.table.ga$False_Negatives)*100
hits.table.ga$Specificity<- (nrow(facs.pcs.top.ga)-hits.table.ga$False_Positives-hits.table.ga$False_Negatives-hits.table.ga$True_Positives)/
  ((nrow(facs.pcs.top.ga)-hits.table.ga$False_Positives-hits.table.ga$False_Negatives-hits.table.ga$True_Positives)+hits.table.ga$False_Positives)*100


write.csv(hits.table.ga, file="ga_EWAS_hit_tables.csv")

