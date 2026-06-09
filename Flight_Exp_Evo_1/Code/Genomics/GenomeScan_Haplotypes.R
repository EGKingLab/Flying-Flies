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

getP<-function(LOD,n,df,ddf)
{
  ff<-(10^((2/n)*LOD)-1)*((ddf)/df)
  pp<- -(pf(ff,df,ddf,lower.tail=FALSE,log=TRUE)/log(10))
  #pp<- pf(ff,df,ddf,lower.tail=FALSE)
  return(pp)
}

getF<-function(LOD,n,df)
{
  ddf<-n-(df+1)
  ff<-(10^((2/n)*LOD)-1)*((n-df-1)/df)
  pp<- -log10(pf(ff,df,ddf,lower.tail=FALSE))
  return(ff)
}

logLik.multi <- function(object, REML = FALSE, ...)
{
  
  all.val<-numeric(length=0)
  for(j in 1:ncol(object$residuals))
  {
    res <- object$residuals[,j] # not resid(object) because of NA methods
    p <- object$rank
    N <- length(res) 
    if(is.null(w <- object$weights)) {
      w <- rep.int(1, N)
    } else {
      ## this is OK as both resids and weights are for the cases used
      excl <- w == 0      # eliminating zero weights
      if (any(excl)) {
        res <- res[!excl]
        N <- length(res)
        w <- w[!excl]
      }
    }
    N0 <- N
    if(REML) N <- N - p
    val <- .5* (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) +
                                     log(sum(w*res^2))))
    if(REML) val <- val - sum(log(abs(diag(object$qr$qr)[1L:p])))
    attr(val, "nall") <- N0 # NB, still omits zero weights
    attr(val, "nobs") <- N
    attr(val, "df") <- p + 1
    class(val) <- "logLik"
    all.val<-rbind(all.val,val)
  }#for close
  all.val
}



library(tidyverse)
library(cowplot)

args <- commandArgs(TRUE)
ChArm <- args[1]

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

output <- data.frame("ChArm" = fdat$intInfo$ChArm,
                     "Ppos" = fdat$intInfo$Ppos,
                     "MA_D" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A1_D" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A2_D" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A3_D" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A4_D" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A5_D" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A6_D" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A7_D" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A8_D" = numeric(length=length(fdat$intInfo$ChArm)),
                     "O_F" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A1_t" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A2_t" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A3_t" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A4_t" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A5_t" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A6_t" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A7_t" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A8_t" = numeric(length=length(fdat$intInfo$ChArm)),
                     "O_p" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A1_p" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A2_p" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A3_p" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A4_p" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A5_p" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A6_p" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A7_p" = numeric(length=length(fdat$intInfo$ChArm)),
                     "A8_p" = numeric(length=length(fdat$intInfo$ChArm)))
                     

haps <- paste0("A",1:8)

for(hapi in haps)
{
 output[,paste0(hapi,"_D")]<- abs(rowMeans(ffreq_all[,hapi,1:4]) - rowMeans(ffreq_all[,hapi,5:8]))

transformed_d <- asin(sqrt(ffreq_all[,hapi,] + 0.001))

group <- rep(c("FP","CON"),each=4)
tt <- rowtstat(X=transformed_d, grp=group)
output[,paste0(hapi, "_t")] <- tt

#adding log(2) does 2 tailed test by doubling p val
output[,paste0(hapi, "_p")] <- -(pt(abs(tt),6,lower.tail=FALSE,log.p=TRUE) + log(2))/log(10)

}

cat("hapscan")

output$MA_D <- rowMeans(output[,paste0(haps,"_D")])

haps.m <- rep(haps, 8)
pool.m <- rep(c("FP","CON"), each = 32)
#t1 <- Sys.time()
mod <- logLik.multi(lm(ffreq_long ~ haps.m + haps.m:pool.m))/log(10)
mod0 <- logLik.multi(lm(ffreq_long ~ haps.m))/log(10)
#t2 <- Sys.time()
LOD <- mod-mod0
output[,"O_p"] <- as.numeric(getP(LOD=LOD, n=64, df = 8, ddf = 48))
#use dspr functions to calculate f from lod
output[,"O_F"] <- as.numeric(getF(LOD=LOD, n=64, df = 8))

cat("overallscan")

#anova(lm(ffreq_long[,100] ~haps.m + haps.m:pool.m))
#close, not exact


saveRDS(output, file = paste0("../../../Big_data/Haplotype_Inf/HapScan_", ChArm, ".Rds"))

