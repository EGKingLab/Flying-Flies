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

Arms <- c("2L","2R","3L","3R","X")
output.all <- data.frame("CHROM"= character(length=0),
                         "POS"= numeric(length=0),
                         "t"= numeric(length=0),
                         "log10P"= numeric(length=0))

for(ii in 1:length(Arms))
{
ss <- readRDS(file = paste0("../../../Big_data/Haplotype_Inf/PoolFreq_EstFrom_Haps_", Arms[ii], ".Rds"))
output <- ss[,1:2]
output$t <- NA
output$log10P <- NA

dd <- ss[,3:10]
dd[dd > 1] <- 1

transformed_d <- asin(sqrt(dd))

group <- rep(c("FP","CON"),each=4)
output[,'t'] <- rowtstat(X=transformed_d, grp=group)

#adding log(2) does 2 tailed test by doubling p val
output[,'log10P'] <- -(pt(abs(output[,'t']),6,lower.tail=FALSE,log.p=TRUE) + log(2))/log(10)
output.all <- rbind(output.all, output)
}

output.all$Pval <- 10^-(output.all$log10P)
output.all$padj <- p.adjust(output.all$Pval, method = "BH")

rawp <- readRDS("/home/srikantv/Flying_Flies_Phenotyping/Big_data/lm_adjpvals_SNPs.RDS")

cc <- left_join(output.all, rawp)

cc$dd <- cc$lm_adjPval - cc$padj

ww <- which(cc$lm_adjPval < 0.05)

cc[ww,]

