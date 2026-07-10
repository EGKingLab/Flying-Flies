library(tidyverse)
library(cowplot)

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

fdat_long <- data.frame("CHROM"= pp[,1],
                        "POS" = as.integer(pp[,2]),
                        "A1"= fdat$mle[,1],
                        "A2"= fdat$mle[,2],
                        "A3"= fdat$mle[,3],
                        "A4"= fdat$mle[,4],
                        "A5"= fdat$mle[,5],
                        "A6"= fdat$mle[,6],
                        "A7"= fdat$mle[,7],
                        "A8"= fdat$mle[,8])








pools <- c("FP1", "FP2", "FP3", "FP4", "CON1","CON2","CON3","CON4")
types <- rep(c("FP","C"),each=4)
  

counter <- 1
for(pool in pools){
  
  fdat <- readRDS(file = paste0("../../../Big_data/Haplotype_Inf/HapInf_", ChArm, "_",
                                pool, ".Rds"))
  fdat_long <- data.frame("CHROM"= fdat$intInfo$ChArm,
                          "POS" = fdat$intInfo$Ppos,
                          "A1"= fdat$mle[,1],
                          "A2"= fdat$mle[,2],
                          "A3"= fdat$mle[,3],
                          "A4"= fdat$mle[,4],
                          "A5"= fdat$mle[,5],
                          "A6"= fdat$mle[,6],
                          "A7"= fdat$mle[,7],
                          "A8"= fdat$mle[,8])
  fdat_long <- pivot_longer(fdat_long, cols= starts_with("A"), names_to = "Haplotype", values_to = "Frequency")
  fdat_long$Pool <- pool
  fdat_long$Type <- types[counter]
  if(counter == 1){
    fdat_all <- fdat_long
  }else{
    fdat_all <- rbind(fdat_all, fdat_long) 
  }
  counter <- counter + 1
}



