rowtstat <- function(X, grp){
  t_stat <- function(X) {
    m <- rowMeans(X)
    n <- ncol(X)
    var <- rowSums((X - m) ^ 2) / (n - 1)
    
    list(m = m, n = n, var = var)
  }
  
  g1 <- t_stat(X[, grp == "FP"])
  g2 <- t_stat(X[, grp == "CON"])
  
  se_total <- sqrt(g1$var / g1$n + g2$var / g2$n)
  (g1$m - g2$m) / se_total
}




library(tidyverse)
library(cowplot)

ChArm <- "2L"
pools <- c("FP1", "FP2", "FP3", "FP4", "CON1","CON2","CON3","CON4")

ffreq_all <- data.frame(
  "sampleID" = character(length = 0),
  "Haplotype" = character(length = 0),
  "Frequency" = numeric(length = 0)
)

fdat <- readRDS(file = paste0("../../../Big_data/Haplotype_Inf/HapInf_", ChArm, "_",
                              pools[1], ".Rds"))
ffreq <- fdat$mle

ffreq_all <- array(, c(nrow(ffreq), ncol(ffreq),length(pools)))
dimnames(ffreq_all)[[1]] <- paste0(fdat$intInfo$ChArm,"_",fdat$intInfo$Ppos)
dimnames(ffreq_all)[[2]] <- paste0("A",1:8)
dimnames(ffreq_all)[[3]] <- pools

ffreq_long <- matrix(,nrow = ncol(ffreq)*length(pools),ncol = nrow(ffreq))
dimnames(ffreq_long)[[1]] <- paste0(rep(paste0("A",1:8),8),"_",rep(pools, each = 8))
dimnames(ffreq_long)[[2]] <- paste0(fdat$intInfo$ChArm,"_",fdat$intInfo$Ppos)
counter <- 1
for(pool in pools){
  
  fdat <- readRDS(file = paste0("../../../Big_data/Haplotype_Inf/HapInf_", ChArm, "_",
                                pool, ".Rds"))
  ffreq_all[,,pool] <- fdat$mle
  
  ffreq_long[counter:(counter+7),] <- t(fdat$mle)
  counter <- counter + 8
  
}

modx <- paste0(rep(paste0("A",1:8),8),"_",rep(c("FP","CON"), each = 32))
t1 <- Sys.time()
mod <- lm(ffreq_long ~ haps.m + haps.m:pool.m)
mod0 <- lm(ffreq_long ~ haps.m)
t2 <- Sys.time()

#use dspr functions to calculate f from lod


temp <- ffreq_all[,1,]


group <- rep(c("FP","CON"),each=4)
tt <- rowtstat(X=temp, grp=group)

#adding log(2) does 2 tailed test by doubling p val
pp <- pt(abs(tt),6,lower.tail=FALSE,log.p=TRUE) + log(2)


