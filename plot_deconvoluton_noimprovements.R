### NOTE this was run locally

load("GenR_cord_deconvolution_predicted_celltypes_noimprovements.rdata")
load("genr_facs_counts.rdata")
rownames(est_cell_counts) <-gsub("/home/040129/3roundEWAs-AllFiles/", "", rownames(est_cell_counts))

validation_counts<- validation_counts[!rownames(validation_counts)%in%c("8667045039_R02C01","7668610116_R03C01", "7668610134_R02C01", "7668610148_R02C02"),]
est_cell_counts<- est_cell_counts[!rownames(est_cell_counts)%in%c("8667045039_R02C01","7668610116_R03C01", "7668610134_R02C01", "7668610148_R02C02"),]
est_cell_counts_orig<- as.data.frame(est_cell_counts)
est_cell_counts<- as.data.frame(est_cell_counts)
est_cell_counts$nRBC<- NULL
dim(est_cell_counts)
dim(validation_counts)

validation_counts<- validation_counts[rownames(est_cell_counts),]


library(reshape2)

est.melt<-melt(cbind(est_cell_counts, ID=rownames(est_cell_counts), type="estimate"), id=c("ID", "type"))
tail(est.melt)
val.melt<- melt(cbind(validation_counts, ID=rownames(validation_counts), type="count"), id=c("ID", "type"))
head(val.melt)

all.melt<- rbind(est.melt, val.melt)
all.cast<- dcast(all.melt, ID+variable~type)
head(all.cast)

ggplot(all.cast, aes(count, estimate))+
  geom_point()+
  facet_wrap(~variable, scales="free")+
  stat_smooth(method="lm")+
  theme_bw()

cor.test(est_cell_counts$Bcell, validation_counts$Bcell, method="spearman")#0.67 (-0.03)
cor.test(est_cell_counts$CD8T, validation_counts$CD8T, method="spearman") #0.68 (+0.08)
cor.test(est_cell_counts$CD4T, validation_counts$CD4T, method="spearman")#0.78 (+0.08)
cor.test(est_cell_counts$Gran, validation_counts$Gran, method="spearman") #0.63 (same)
cor.test(est_cell_counts$Mono, validation_counts$Mono, method="spearman") #0.4 (+0.02)
cor.test(est_cell_counts$NK, validation_counts$NK, method="spearman") #0.8 (same)

mad(est_cell_counts$Gran) #0.09
mad(est_cell_counts$CD4T) #0.04
mad(est_cell_counts$CD8T) #0.017
mad(est_cell_counts$Bcell) #0.014
mad(est_cell_counts$NK) #0.03
mad(est_cell_counts$Mono) #0.015

ggplot(est_cell_counts_orig, aes(est_cell_counts_orig$nRBC))+
  geom_histogram()+
  theme_bw()

cor.mad.noimp<-  data.frame(matrix(NA, 7, 2, dimnames=list(c("CD8T","CD4T","NK","Bcell","Mono", "Gran", "nRBC")
                                                 ,c("Rho", "MAD"))),stringsAsFactors=F)

for (i in c("CD8T","CD4T","NK","Bcell","Mono", "Gran", "nRBC")){
  cor.mad.noimp[i,"MAD"]<- mad(est_cell_counts_orig[,i])
  cor.mad.noimp[i,"Rho"]<- cor.test(est_cell_counts[,i], validation_counts[,i], method="spearman")$estimate
}

save(cor.mad.noimp, file="noimp_predictions_cor_and_mad.rdata")
