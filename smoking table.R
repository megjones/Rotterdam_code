setwd("C:/Users/040129/Documents/Genr/toptables")

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
hits.table.smoke<- data.frame(matrix(NA, 10, 3,
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
  hits.table.smoke[a, 1]<- length(intersect(rownames(top)[top$adj.P.Val< 0.2], 
                                            rownames(top)[abs(top$db)>0.02]))
  hits.table.smoke[a, 2]<- length(intersect(intersect(rownames(top)[top$adj.P.Val< 0.2], 
                                                      rownames(top)[abs(top$db)>0.02]),true_hits_smoke))
  hits.table.smoke[a,3]<- hits.table.smoke[a,1]-hits.table.smoke[a,2]
}

write.csv(hits.table.smoke, file="smoke_EWAS_hit_tables.csv")

