


ChArms <- c("X","2L","2R","3L","3R")
pools <- c("FP1", "FP2", "FP3", "FP4", "CON1","CON2","CON3","CON4")

for(ChArm in ChArms){
fdat <- readRDS(file = paste0("../../../Big_data/Haplotype_Inf/HapInf_", ChArm, "_",
                              pools[1], ".Rds"))
ffreq <- fdat$mle

ffreq_all <- array(, c(nrow(ffreq), ncol(ffreq),length(pools)))
dimnames(ffreq_all)[[1]] <- paste0(fdat$intInfo$ChArm,"_",fdat$intInfo$Ppos)
dimnames(ffreq_all)[[2]] <- paste0("A",1:8)
dimnames(ffreq_all)[[3]] <- pools

intInfo <- fdat$intInfo

for(pool in pools){
  
  fdat <- readRDS(file = paste0("../../../Big_data/Haplotype_Inf/HapInf_", ChArm, "_",
                                pool, ".Rds"))
  ffreq_all[,,pool] <- fdat$mle
  
}

intInfo$FP1_sum <- rowSums(ffreq_all[,,1])
intInfo$FP2_sum <- rowSums(ffreq_all[,,2])
intInfo$FP3_sum <- rowSums(ffreq_all[,,3])
intInfo$FP4_sum <- rowSums(ffreq_all[,,4])
intInfo$CON1_sum <- rowSums(ffreq_all[,,5])
intInfo$CON2_sum <- rowSums(ffreq_all[,,6])
intInfo$CON3_sum <- rowSums(ffreq_all[,,7])
intInfo$CON4_sum <- rowSums(ffreq_all[,,8])

saveRDS(intInfo, file = paste0("../../../Big_data/Haplotype_Inf/IntInf_", ChArm, ".Rds"))

saveRDS(ffreq_all, file = paste0("../../../Big_data/Haplotype_Inf/Array_HapFreq_", ChArm, ".Rds"))
}
