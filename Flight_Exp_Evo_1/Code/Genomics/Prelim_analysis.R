
ss <- read.table("../../../Big_data/FlyingFlies_SNPtable.txt", header=TRUE)

#grab frequencies only - in every other column starting with 6
#control are 1st 4, fp are second
ss_freq <- ss[,seq(6,20,by=2)]

ss_freq <- as.matrix(ss_freq)

#arcsin sqrt transform
ss_freq <- asin(sqrt(ss_freq))

output <- numeric(length = nrow(ss_freq))


for(ii in 1:nrow(ss_freq))
{


ss.i <- data.frame("Treatment" = rep(c("C","FP"), each = 4), 
                   "Frequency" = ss_freq[ii,])

mod0 <- lm(Frequency ~ 1, data = ss.i)
mod1 <- lm(Frequency ~ Treatment, data = ss.i)

aa <- anova(mod0, mod1)

output[ii] <- aa$`Pr(>F)`[2]

cat(ii,"\n")
}


final.out <- ss[,c(2,3)]
final.out$Pval <- output
final.out$Pval_adj <- p.adjust(output, method = "fdr")

saveRDS(final.out, file = "../../../Big_data/Prelim_Scan.RDS")

rr <- readRDS("../../../Big_data/Prelim_Scan.RDS")
