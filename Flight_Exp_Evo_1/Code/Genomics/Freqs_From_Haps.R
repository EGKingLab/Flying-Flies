
library(tidyverse)
library(readxl)

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


## Base pop


rr <- read_excel("../../../Raw_data/RIL_IDs_FF_Hamlin_etal_2025.xlsx")|>
  mutate(ril = as.integer(ril), rilC = as.character(ril))

rilsnp <- readRDS(file="/home/kingeg/Projects/MyGithub/HaplotypeEstimator/BigData/SNPtable_A_RILs_inferred_Release4Conv6_X.rds")
rilsnp <- rilsnp[,rr$rilC]
pp <- str_split(rownames(rilsnp),"_", simplify = TRUE)

afreq_all <- tibble(CHROM = pp[,1],
                POS = as.integer(pp[,2]),
                Base_freq = rowMeans(rilsnp))

for(Arm in c("2L","2R","3L","3R"))
{
rilsnp <- readRDS(file=paste0("/home/kingeg/Projects/MyGithub/HaplotypeEstimator/BigData/SNPtable_A_RILs_inferred_Release4Conv6_",Arm,".rds"))
rilsnp <- rilsnp[,rr$rilC]
pp <- str_split(rownames(rilsnp),"_", simplify = TRUE)

afreq <- tibble(CHROM = pp[,1],
                POS = as.integer(pp[,2]),
                Base_freq = rowMeans(rilsnp))
afreq_all <- rbind(afreq_all, afreq)
}


filt.snp <- readRDS(file = "/home/srikantv/Flying_Flies_Phenotyping/Big_data/FF_genomes_filtered.RDS")

base_freq <- left_join(filt.snp, afreq_all)

#write out just base
saveRDS(base_freq[,c("CHROM","POS","Base_freq")], file = "../../../Big_data/Base_pop_SNP_freq.RDS")
