library(ggridges)

myColors <- c("gold","goldenrod",
              "#2171b5","#6baed6","#6baed6",
              "#9467BDFF","#D62728FF",
              "#238443","#78c679","#addd8e","#d9f0a3",
              "#dd3497","#fa9fb5","grey")

color_possibilities<-c("FACS - PCA - Gold-Standard","FACS - Drop One Cell Type",
                       "Deconvolution - PCA","Deconvolution - Drop One Cell Type","Deconvolution - Counts",
                       "ReFACTor","RefFreeCellMix",
                       "SVA - Supervised Smoking","SVA - Unsupervised Smoking",
                       "SVA - Supervised Sex","SVA - Unsupervised Sex",
                       "RUV - Smoking", "RUV - Sex",
                       "Uncorrected")

names(myColors) <- color_possibilities
fillscale <- scale_fill_manual(name="Method",
                               values = myColors, drop = T)

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

smoke.p<- data.frame(pval=c(facs.pcs.top.smoke$P.Value,
                            facs.top.smoke$P.Value,
                            decon.pcs.top.smoke$P.Value,
                            decon.top.smoke$P.Value,
                            refactor.top.smoke$P.Value,
                            reffreecellmix.top.smoke$P.Value,
                            sva.sup.top.smoke$P.Value,
                            sva.unsup.top.smoke$P.Value,
                            ruv.top.smoke$P.Value,
                            uncorrected.top.smoke$P.Value
                            ),
                     Method=c(rep("FACS - PCA - Gold-Standard", nrow(facs.pcs.top.smoke)),
                              rep("FACS - Drop One Cell Type", nrow(facs.top.smoke)),
                              rep("Deconvolution - PCA", nrow(decon.pcs.top.smoke)),
                              rep("Deconvolution - Drop One Cell Type", nrow(decon.top.smoke)),
                              rep("ReFACTor", nrow(refactor.top.smoke)),
                              rep("RefFreeCellMix", nrow(reffreecellmix.top.smoke)),
                              rep("SVA - Supervised Smoking", nrow(sva.sup.top.smoke)),
                              rep("SVA - Unsupervised Smoking", nrow(sva.unsup.top.smoke)),
                              rep("RUV - Smoking", nrow(ruv.top.smoke)),
                              rep("Uncorrected", nrow(uncorrected.top.smoke))
                              )) 
ggplot(smoke.p, aes(pval,  fill=Method)) +
  geom_density(aes(alpha=0.5)) + 
  fillscale   +xlab("Unadjusted Pvalue distribution")+ylab("")


smoke.plot<- ggplot(smoke.p, aes(pval,  fill=Method)) +
  geom_density(aes(alpha=0.5)) + 
  fillscale   +xlab("Unadjusted Pvalue distribution")+ylab("")+
  theme(legend.position="none")

### repeat for sex

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

sex.p<- data.frame(pval=c(facs.pcs.top.sex$P.Value,
                            facs.top.sex$P.Value,
                            decon.pcs.top.sex$P.Value,
                            decon.top.sex$P.Value,
                            refactor.top.sex$P.Value,
                            reffreecellmix.top.sex$P.Value,
                            sva.sup.top.sex$P.Value,
                            sva.unsup.top.sex$P.Value,
                            ruv.top.sex$P.Value,
                            uncorrected.top.sex$P.Value
),
Method=c(rep("FACS - PCA - Gold-Standard", nrow(facs.pcs.top.sex)),
         rep("FACS - Drop One Cell Type", nrow(facs.top.sex)),
         rep("Deconvolution - PCA", nrow(decon.pcs.top.sex)),
         rep("Deconvolution - Drop One Cell Type", nrow(decon.top.sex)),
         rep("ReFACTor", nrow(refactor.top.sex)),
         rep("RefFreeCellMix", nrow(reffreecellmix.top.sex)),
         rep("SVA - Supervised GA", nrow(sva.sup.top.sex)),
         rep("SVA - Unsupervised GA", nrow(sva.unsup.top.sex)),
         rep("RUV - GA", nrow(ruv.top.sex)),
         rep("Uncorrected", nrow(uncorrected.top.sex))
)) 

ggplot(sex.p, aes(pval,  fill=Method)) +
  geom_density(aes(alpha=0.5)) + 
  fillscale   +xlab("Unadjusted Pvalue distribution")+ylab("")

sex.plot<- ggplot(sex.p, aes(pval,  fill=Method)) +
  geom_density(aes(alpha=0.5)) + 
  fillscale   +xlab("Unadjusted Pvalue distribution")+ylab("")+
  theme(legend.position="none")

## repeat for GA

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

ga.p<- data.frame(pval=c(facs.pcs.top.ga$P.Value,
                          facs.top.ga$P.Value,
                          decon.pcs.top.ga$P.Value,
                          decon.top.ga$P.Value,
                          refactor.top.ga$P.Value,
                          reffreecellmix.top.ga$P.Value,
                          sva.sup.top.ga$P.Value,
                          sva.unsup.top.ga$P.Value,
                          ruv.top.ga$P.Value,
                          uncorrected.top.ga$P.Value
),
Method=c(rep("FACS - PCA - Gold-Standard", nrow(facs.pcs.top.ga)),
         rep("FACS - Drop One Cell Type", nrow(facs.top.ga)),
         rep("Deconvolution - PCA", nrow(decon.pcs.top.ga)),
         rep("Deconvolution - Drop One Cell Type", nrow(decon.top.ga)),
         rep("ReFACTor", nrow(refactor.top.ga)),
         rep("RefFreeCellMix", nrow(reffreecellmix.top.ga)),
         rep("SVA - Supervised GA", nrow(sva.sup.top.ga)),
         rep("SVA - Unsupervised GA", nrow(sva.unsup.top.ga)),
         rep("RUV - GA", nrow(ruv.top.ga)),
         rep("Uncorrected", nrow(uncorrected.top.ga))
)) 



ga.plot<- ggplot(ga.p, aes(pval,  fill=Method)) +
  geom_density(aes(alpha=0.5)) + 
  fillscale   +xlab("Unadjusted Pvalue distribution")+ylab("")+
  theme(legend.position="none")


grid.arrange(sex.plot, ga.plot, smoke.plot, nrow=1)
