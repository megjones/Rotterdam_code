mads.melt<- melt(mads, id="celltypes")
ggplot(mads.melt, aes(celltypes, variable, fill=value))+
  geom_tile()+
  geom_text(aes(label = round(value, 3)), size=3) +
  theme_bw()+
  scale_fill_gradient(low = "white",  high = "cornflowerblue")

load("~/Code for ref free methods/Louie_predicted_WB_celltypes_Oct27_use_Louies_sorted.rdata")
uncorrected<- est_cell_counts
plot.counts<- melt(rbind(cbind(as.data.frame(uncorrected), Method="uncorrected"),
                         cbind(as.data.frame(facs.counts), Method="FACS_counts"), 
                         cbind(as.data.frame(facs.pcs.counts), Method="FACS_PCs"), 
                         cbind(as.data.frame(decon.counts), Method="Decon_counts"), 
                         cbind(as.data.frame(decon.pcs.counts), Method="Decon_PCs"), 
                         cbind(as.data.frame(sva.unsup.sex), Method="sva.unsup.sex"), 
                         cbind(as.data.frame(sva.unsup.ga), Method="sva.unsup.ga"), 
                         cbind(as.data.frame(sva.sup.ga), Method="sva.sup.ga"), 
                         cbind(as.data.frame(ruv.sex), Method="ruv.sex"), 
                         cbind(as.data.frame(ruv.ga), Method="ruv.ga"), 
                         cbind(as.data.frame(cellmix), Method="cellmix"), 
                         cbind(as.data.frame(refactor), Method="refactor")
), id="Method")


#ggplot(plot.counts, aes(Method, value))+
#  geom_density_ridges()+
#  facet_wrap(~variable)+
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#ggplot(plot.counts, aes(Method, value))+
#  geom_violin()+
#  facet_wrap(~variable)+
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(plot.counts, aes(Method, value))+
  geom_point(position = position_jitter(width = 0.1))+
  facet_wrap(~variable)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


