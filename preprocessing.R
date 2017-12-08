# Ripped from ChAMP
preprocess_rgset <- function(rgSet, detPcut = 0.01, beadCutoff=0.05, methValue="B", 
                             preprocess_method=preprocessRaw) {
  sampleNames(rgSet) = rgSet[[1]]
  pd <- pData(rgSet)
  green = getGreen(rgSet)
  red = getRed(rgSet)
  mset <- preprocess_method(rgSet)
  detP <- detectionP(rgSet)
  failed <- detP > detPcut
  
  # Failed detection probes.
  numfail = colMeans(failed)
  message("The fraction of failed positions per sample: ")
  print(numfail)
  
  # Filter failed probes.
  # Fraction of allowed failed probes.
  removeDetP = 0
  mset.f = mset[rowSums(detP >= detPcut) <= removeDetP *  ncol(detP), ]
  if (removeDetP == 0) {
    message("Filtering probes with a detection p-value above ", 
            detPcut, " in one or more samples has removed ", 
            dim(mset)[1] - dim(mset.f)[1], " probes from the analysis. If a large number of probes have been removed, ChAMP suggests you look at the failedSample.txt file to identify potentially bad samples.")
  } else {
    message("Filtering probes with a detection p-value above ", 
            detPcut, " in at least ", removeDetP * 100, "% of samples has removed ", 
            dim(mset)[1] - dim(mset.f)[1], " probes from the analysis. If a large number of probes have been removed, ChAMP suggests you look at the failedSample.txt file to identify potentially bad samples.")
  }
  mset = mset.f
  
  # Filter for low bead coverage if available.

    bc = beadcount(rgSet)
    bc2 = bc[rowSums(is.na(bc)) < beadCutoff * (ncol(bc)), ]
    mset.f2 = mset[featureNames(mset) %in% row.names(bc2), ]
    message("Filtering probes with a beadcount <3 in at least ", 
            beadCutoff * 100, "% of samples, has removed ", dim(mset)[1] - 
              dim(mset.f2)[1], " from the analysis.")
    mset = mset.f2
  
  
  # Filter SNPs
  data(snp.hit)
  mset.f2 = mset[!featureNames(mset) %in% snp.hit$TargetID,  ]
  message("Filtering probes with SNPs as identified in Nordlund et al, has removed ", 
          dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
  mset = mset.f2
  
  # Filter cross-hybridizing probes.
  data(multi.hit)
  mset.f2 = mset[!featureNames(mset) %in% multi.hit$TargetID,  ]
  message("Filtering probes that align to multiple locations as identified in Nordlund et al, has removed ", 
          dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
  mset = mset.f2
  
  # Get Betas.
  intensity = minfi::getMeth(mset) + minfi::getUnmeth(mset)
  beta.raw = minfi::getBeta(mset, "Illumina")
  detP = detP[which(row.names(detP) %in% row.names(beta.raw)),  ]
  
  # Filter out Sex Chromosome Probes.
  data(probe.features)
  autosomes = probe.features[!probe.features$CHR %in% c("X","Y"), ]
  mset.f2 = mset[featureNames(mset) %in% row.names(autosomes),  ]
  beta.raw = beta.raw[row.names(beta.raw) %in% row.names(autosomes),  ]
  detP = detP[row.names(detP) %in% row.names(autosomes),  ]
  intensity = intensity[row.names(intensity) %in% row.names(autosomes),  ]
  message("Filtering probes on the X or Y chromosome has removed ", 
          dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
  mset = mset.f2
  
  message("The analysis will proceed with ", dim(beta.raw)[1], 
          " probes and ", dim(beta.raw)[2], " samples.")
  
  mset
}