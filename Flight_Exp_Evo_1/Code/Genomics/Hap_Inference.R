library(tidyverse)
library(INLA)

source("../../../Functions/Hap_functions.R")
source("../../../Functions/read_from_h5.R")


#ChArm <- "2L"
#pool <- "FP1" #column of data

args <- commandArgs(TRUE)
ChArm <- args[1]
pool <- args[2]

run.type <- "ppos"
Nhaps <- 8

pool.snp <- readRDS(file= "../../../Big_data/SNPpool_wFounders.RDS")

pool.snp <- pool.snp |>
  filter(CHROM == ChArm)

p.names <- paste0(pool.snp$CHROM,"_", pool.snp$POS)

Foundsnp <- select(pool.snp, contains("freq_A"))
rownames(Foundsnp) <- p.names

pool.snp <- select(pool.snp, contains(pool))
pool.snp <- as.matrix(pool.snp)
rownames(pool.snp) <- p.names

poslist <- readRDS(paste0("../../../Big_data/Pos_Phys_Gen_",ChArm,".rds"))

ww <- which(poslist$POSCODE %in% rownames(pool.snp))
poslist <- poslist[ww,]

pool.snp <- pool.snp[poslist$POSCODE,]
Foundsnp <- Foundsnp[poslist$POSCODE,]




test.positions <- 1:nrow(poslist)
CC <- pool.snp[,2]
newf <- pool.snp[,1]

#set up output holding
counter <- 1
mle.res <- matrix(NA, nrow = length(test.positions), ncol = Nhaps)
#inla.res <- matrix(NA, nrow = length(test.positions), ncol = Nhaps)

interval_info_all <- data.frame(
  "ChArm" = character(length=length(test.positions)),
  "Ppos" = numeric(length=length(test.positions)),
  "type" = character(length=length(test.positions)),
  "nsnp" = numeric(length=length(test.positions)),
  "Pdist" = numeric(length=length(test.positions)),
  "Gdist" = numeric(length=length(test.positions)),
  "cormax" = numeric(length=length(test.positions)),
  "nHcor" = numeric(length=length(test.positions))
)


initials.mle <- as.matrix(rep(1/Nhaps, Nhaps), ncol = 1)
inla_prior_mean <- list(default = 0.125)
# Set prior for INLA precision
inla_prior_precision <- list(default = precision(0.5))

for(ii in test.positions){
  
  dat <- get_pos(ii, newf, Foundsnp, Nhaps,
                 100000, CC = CC, type="ppos", poslist)
  
  
  # MLE
  mle.hap <- est.mle(X = dat$X, Y = dat$Y,
                     p = initials.mle, CC = dat$Ccov)
  
  # INLA
  
  # fm_inla <- inla(
  #   Y ~ -1 + freq_A1 + freq_A2 + freq_A3 + freq_A4 +
  #     freq_A5 + freq_A6 + freq_A7 + freq_AB8,
  #   control.fixed = list(mean = inla_prior_mean,
  #                        prec = inla_prior_precision),
  #   data = data.frame(Y = dat$Y, dat$X),
  #   num.threads = 1)
  # 
  # inla.hap <- fm_inla$summary.fixed$mean
  
  # Assign output
  mle.res[counter,] <- mle.hap
  #inla.res[counter,] <- inla.hap
  interval_info_all[counter,] <- dat$IntInfo
  counter <- counter + 1
  cat(counter, "\n")
}


ll <- list("intInfo"= interval_info_all,
           "mle" = mle.res)
           #"inla" = inla.res)

saveRDS(ll, file = paste0("../../../Big_data/Haplotype_Inf/HapInf_", ChArm, "_",
                          pool, ".Rds"))
      
