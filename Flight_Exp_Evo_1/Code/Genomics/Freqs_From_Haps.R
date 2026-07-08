
library(tidyverse)

pool.snp <- readRDS(file= "../../../Big_data/SNPpool_wFounders.RDS")

Arms <- c("X","2L","2R","3L","3R")
pools <- c("FP1", "FP2", "FP3", "FP4", "CON1","CON2","CON3","CON4")

for(ii in 1:length(Arms))
{
  
fdat <- readRDS(file = paste0("../../../Big_data/Haplotype_Inf/HapInf_", Arms[ii], "_",
                              pools[1], ".Rds"))

pos.list <- fdat$intInfo[,c("ChArm","Ppos")]
colnames(pos.list) <- c("CHROM","POS")
hdat <- inner_join(pool.snp, pos.list)

est_pool_all <- matrix(data = NA,nrow = nrow(pos.list), ncol = length(pools))

for(jj in 1:length(pools))
{
  fdat <- readRDS(file = paste0("../../../Big_data/Haplotype_Inf/HapInf_", Arms[ii], "_",
                                pools[jj], ".Rds"))
  ffreq <- fdat$mle
  est_pool_all[,jj] <- rowSums(ffreq*hdat[,c("freq_A1","freq_A2","freq_A3","freq_A4",
                                    "freq_A5","freq_A6","freq_A7","freq_AB8")])
  
}
colnames(est_pool_all) <- pools
est_pool_all <- cbind(pos.list, est_pool_all)
saveRDS(est_pool_all, file = paste0("../../../Big_data/Haplotype_Inf/PoolFreq_EstFrom_Haps_", Arms[ii], ".Rds"))
}




