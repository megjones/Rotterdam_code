
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

