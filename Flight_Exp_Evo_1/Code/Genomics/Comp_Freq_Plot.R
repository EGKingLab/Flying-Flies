library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())


ChArm <- "2L"

#get mean freq per haplotype
#get difference from mean for each 
#put in long format

ff <- readRDS(file = paste0("../../../Big_data/Haplotype_Inf/Array_HapFreq_", ChArm, ".Rds"))
ffdiff <- ff

for(ii in 1:8){

pmeans <- rowMeans(ff[,ii,])
ffdiff[,ii,] <- ff[,ii,]-pmeans

}

pp <- str_split(dimnames(ffdiff)[[1]], "_", simplify = TRUE)

pp <- pp[10000:105000,]
ffdiff <- ffdiff[10000:105000,,]

fdat_long <- data.frame("CHROM"= rep(pp[,1],times = 64),
                        "POS" = rep(as.integer(pp[,2]),times = 64),
                        "Pool" = rep(dimnames(ffdiff)[[3]],each=nrow(pp), times= 8),
                        "Linetype" = rep(c("FP","CON"),each = nrow(pp)*4,times=8),
                        "Haplotype" = rep(dimnames(ffdiff)[[2]], each = 8*nrow(pp)),
                        "Freq" = c(ffdiff[,1,1],
                                   ffdiff[,1,2],
                                   ffdiff[,1,3],
                                   ffdiff[,1,4],
                                   ffdiff[,1,5],
                                   ffdiff[,1,6],
                                   ffdiff[,1,7],
                                   ffdiff[,1,8],
                                   ffdiff[,2,1],
                                   ffdiff[,2,2],
                                   ffdiff[,2,3],
                                   ffdiff[,2,4],
                                   ffdiff[,2,5],
                                   ffdiff[,2,6],
                                   ffdiff[,2,7],
                                   ffdiff[,2,8],
                                   ffdiff[,3,1],
                                   ffdiff[,3,2],
                                   ffdiff[,3,3],
                                   ffdiff[,3,4],
                                   ffdiff[,3,5],
                                   ffdiff[,3,6],
                                   ffdiff[,3,7],
                                   ffdiff[,3,8],
                                   ffdiff[,4,1],
                                   ffdiff[,4,2],
                                   ffdiff[,4,3],
                                   ffdiff[,4,4],
                                   ffdiff[,4,5],
                                   ffdiff[,4,6],
                                   ffdiff[,4,7],
                                   ffdiff[,4,8],
                                   ffdiff[,5,1],
                                   ffdiff[,5,2],
                                   ffdiff[,5,3],
                                   ffdiff[,5,4],
                                   ffdiff[,5,5],
                                   ffdiff[,5,6],
                                   ffdiff[,5,7],
                                   ffdiff[,5,8],
                                   ffdiff[,6,1],
                                   ffdiff[,6,2],
                                   ffdiff[,6,3],
                                   ffdiff[,6,4],
                                   ffdiff[,6,5],
                                   ffdiff[,6,6],
                                   ffdiff[,6,7],
                                   ffdiff[,6,8],
                                   ffdiff[,7,1],
                                   ffdiff[,7,2],
                                   ffdiff[,7,3],
                                   ffdiff[,7,4],
                                   ffdiff[,7,5],
                                   ffdiff[,7,6],
                                   ffdiff[,7,7],
                                   ffdiff[,7,8],
                                   ffdiff[,8,1],
                                   ffdiff[,8,2],
                                   ffdiff[,8,3],
                                   ffdiff[,8,4],
                                   ffdiff[,8,5],
                                   ffdiff[,8,6],
                                   ffdiff[,8,7],
                                   ffdiff[,8,8]
                                   ))


hplot <- ggplot(fdat_long, aes(POS,Freq, color = Linetype)) +
  geom_point(size = 0.2, alpha = 1/4) +
  facet_grid(Haplotype ~.)

#add snp
#add sig spots. linetype for hap vs snp

ggsave(hplot, file = paste0("../../Plots/Comparison_Scan_",ChArm,".png"), width = 6.5, height = 6.5)




