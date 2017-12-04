## Preprocessing code - noob/BMIQ
library(minfi)
library(sva)

# Needed for bad probe list.
library(ChAMP)
library(wateRmelon)
library(limma)
library(IlluminaHumanMethylation450kmanifest)
setwd("/home/040129/UBC_files")

genr_sheet <- read.csv("~/3roundEWAs-AllFiles/samplesSheet.csv")
genr_sheet$Basename<-sapply(1:nrow(genr_sheet), function(x) paste("/home/040129/3roundEWAs-AllFiles/", genr_sheet$Basename[x], sep=""))
genr_rgset <- read.metharray.exp(targets=genr_sheet, extended=TRUE) 
source("~/UBC_files/preprocessing.R")
genr_mset.noob <- preprocess_rgset(genr_rgset, preprocess_method = preprocessNoob)
genr_betas.noob <- getBeta(genr_mset.noob)
genr_pd.noob <- pData(genr_mset.noob)


type1<- as.data.frame(getProbeInfo(IlluminaHumanMethylation450kmanifest, type="I"))
type2<- as.data.frame(getProbeInfo(IlluminaHumanMethylation450kmanifest, type="II"))
type<- data.frame(CpG=c(type1$Name, type2$Name), type=c(rep("1", nrow(type1)), rep("2", nrow(type2))))
rownames(type)<- type$CpG
type<- type[rownames(genr_betas.noob),]
type<- as.numeric(type$type)

genr.noob.bmiq<- apply(genr_betas.noob, 2, function(x) BMIQ(x, type)) 

## have to extract the information from this weird format
genr.betas.noob.bmiq<-matrix(NA, nrow=nrow(genr_betas.noob), ncol=ncol(genr_betas.noob))
rownames(genr.betas.noob.bmiq)<- rownames(genr_betas.noob)
colnames(genr.betas.noob.bmiq)<- colnames(genr_betas.noob)
for(i in 1:ncol(genr.betas.noob.bmiq)){
  genr.betas.noob.bmiq[,i]<- genr.noob.bmiq[[i]]$nbeta
}


save(genr.betas.noob.bmiq, file="genr_betas_noob_BMIQ.rdata")
save(genr_pd.noob, file="genr_pdata.rdata")
