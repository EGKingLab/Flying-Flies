library(tidyverse)
library(ggrastr)
library(cowplot)
theme_set(theme_cowplot())

Reg_plot <- function(ChArm, StPos, EnPos)
{
  psnp <- readRDS(paste0("/home/srikantv/Flying_Flies_Phenotyping/Big_data/lm_targetSNPs.RDS"))
  
  psnp <- psnp |>
    filter(CHROM == ChArm & POS >= StPos & POS <= EnPos)
  
  #get mean freq per haplotype
  #get difference from mean for each 
  #put in long format
  
  ff <- readRDS(file = paste0("../../../Big_data/Haplotype_Inf/Array_HapFreq_", ChArm, ".Rds"))
  
  
  posL <- as.integer(str_split(dimnames(ff)[[1]], "_", simplify = TRUE)[,2])
  posL <- posL[posL >= StPos & posL <= EnPos]
  ff <- ff[paste0(ChArm,"_",posL),,]
  
  ffdiff <- ff
  
  for(ii in 1:8){
    
    pmeans <- rowMeans(ff[,ii,])
    ffdiff[,ii,] <- ff[,ii,]-pmeans
    
  }
  
  pp <- str_split(dimnames(ffdiff)[[1]], "_", simplify = TRUE)
  
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
  
  
  
  #add snp
  
  pool.snp <- readRDS(file = "/home/srikantv/Flying_Flies_Phenotyping/Big_data/FF_genomes_filtered.RDS")
  pool.snp <- pool.snp |>
    filter(CHROM == ChArm & POS >= StPos & POS <= EnPos)
  
  freq.snp <- pool.snp |>
    dplyr::select(starts_with("freq_"))
  
  mm <- rowMeans(freq.snp)
  
  #checked this does the diff correct
  sc.freq <- freq.snp - mm
  
  snp_rb <- data.frame("CHROM"= rep(pool.snp$CHROM, times = 8),
                       "POS" = rep(pool.snp$POS, times = 8),
                       "Pool" = rep(dimnames(ffdiff)[[3]],each=nrow(pool.snp)),
                       "Linetype" = rep(c("FP","CON"),each = nrow(pool.snp)*4),
                       "Haplotype" = rep("SNP",nrow(pool.snp)*8),
                       "Freq" = c(sc.freq[,5],
                                  sc.freq[,6],
                                  sc.freq[,7],
                                  sc.freq[,8],
                                  sc.freq[,1],
                                  sc.freq[,2],
                                  sc.freq[,3],
                                  sc.freq[,4]))
  
  fdat_long <- rbind(fdat_long, snp_rb)
  
  phap <- readRDS(file = paste0("../../../Big_data/Haplotype_Inf/HapScan_Padj.Rds"))
  phap <- phap[phap$ChArm==ChArm & phap$Ppos >= StPos & phap$Ppos <= EnPos & phap$Padj <= 0.05,]
  
  
  hplot <- ggplot(fdat_long, aes(POS/1e6,Freq, color = Linetype)) +
    #geom_vline(xintercept = phap$Ppos/1e6, linewidth = 0.5, color = "grey20") +
    #geom_vline(xintercept = psnp$POS/1e6, linewidth = 0.5, color = "gold") +
    geom_point(size = 0.5, alpha = 1/6) +
    facet_grid(Haplotype ~.) +
    xlab("Position (MB)") +
    ylab("Frequency") +
    scale_color_manual(name = "Selection",
                       values = c("CON" = 'steelblue', "FP" = 'firebrick'))
  return(hplot)
}


ChArm <- "3R"
#make sure to use as.integer so it prints in not sci notation
StPos <- as.integer(10500000)
EnPos <- as.integer(11500000)

p1 <- Reg_plot(ChArm, StPos, EnPos)

#I recommend saving to file, not trying to look at it in the window- too big
ggsave(p1, file = paste0("../../Plots/Comparison_Scan_",ChArm,"_",StPos,"_",EnPos,".png"), width = 6.5, height = 6.5)

ChArm <- "3R"
#make sure to use as.integer so it prints in not sci notation
StPos <- as.integer(10000000)
EnPos <- as.integer(15000000)

p1 <- Reg_plot(ChArm, StPos, EnPos)

#I recommend saving to file, not trying to look at it in the window- too big
ggsave(p1, file = paste0("../../Plots/Comparison_Scan_",ChArm,"_",StPos,"_",EnPos,".png"), width = 6.5, height = 6.5)


ChArm <- "3R"
#make sure to use as.integer so it prints in not sci notation
StPos <- as.integer(27969000)
EnPos <- as.integer(27975000)

p1 <- Reg_plot(ChArm, StPos, EnPos)

#I recommend saving to file, not trying to look at it in the window- too big
ggsave(p1, file = paste0("../../Plots/Comparison_Scan_",ChArm,"_",StPos,"_",EnPos,".png"), width = 4.5, height = 5.5)

ChArm <- "3R"
#make sure to use as.integer so it prints in not sci notation
StPos <- as.integer(11011656)
EnPos <- as.integer(11015170)

p1 <- Reg_plot(ChArm, StPos, EnPos)

#I recommend saving to file, not trying to look at it in the window- too big
ggsave(p1, file = paste0("../../Plots/Comparison_Scan_",ChArm,"_",StPos,"_",EnPos,".png"), width = 4.5, height = 5.5)

