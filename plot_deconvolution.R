### NOTE this was run locally

setwd("C:/Users/040129/Documents/Genr")

#validation_counts<- read.table("AnalysisfileUBC.dat", sep="\t", header=T)
#metadata<- validation_counts[!rownames(validation_counts)%in%c("8667045039_R02C01","7668610116_R03C01", "7668610134_R02C01", "7668610148_R02C02"),1:2]
#colnames(metadata)<- c("Sex", "GA")
#metadata$Sex<- as.factor(metadata$Sex)
#save(metadata, file="genr_meta.rdata")

# rownames(validation_counts)<- validation_counts$Sample_ID
# validation_counts<- validation_counts[,7:12]
# head(validation_counts)
# validation_counts<- validation_counts[rownames(est_cell_counts),]
# colnames(validation_counts)<- c("NK", "Bcell", "CD8T", "CD4T", "Gran", "Mono")
# dim(validation_counts)
# validation_counts<- validation_counts[!rownames(validation_counts)%in%c("8667045039_R02C01","7668610116_R03C01", "7668610134_R02C01", "7668610148_R02C02"),]
# save(validation_counts, file="genr_facs_counts.rdata")
load("C:/Users/040129/Documents/Genr/GenR_cord_deconvolution_predicted_celltypes.rdata")
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

cor.test(est_cell_counts$Bcell, validation_counts$Bcell, method="spearman")#0.71
cor.test(est_cell_counts$CD8T, validation_counts$CD8T, method="spearman") #0.62
cor.test(est_cell_counts$CD4T, validation_counts$CD4T, method="spearman")#0.75
cor.test(est_cell_counts$Gran, validation_counts$Gran, method="spearman") #0.63 
cor.test(est_cell_counts$Mono, validation_counts$Mono, method="spearman") #0.42
cor.test(est_cell_counts$NK, validation_counts$NK, method="spearman") #0.81


ggplot(est_cell_counts_orig, aes(est_cell_counts_orig$nRBC))+
  geom_histogram()+
  theme_bw()

cor.mad<-  data.frame(matrix(NA, 7, 2, dimnames=list(c("CD8T","CD4T","NK","Bcell","Mono", "Gran", "nRBC")
                                                           ,c("Rho", "MAD"))),stringsAsFactors=F)

for (i in c("CD8T","CD4T","NK","Bcell","Mono", "Gran", "nRBC")){
  cor.mad[i,"MAD"]<- mad(est_cell_counts_orig[,i])
  cor.mad[i,"Rho"]<- cor.test(est_cell_counts[,i], validation_counts[,i], method="spearman")$estimate
}

save(cor.mad, file="predictions_cor_and_mad.rdata")

difs<- cor.mad-cor.mad.noimp

