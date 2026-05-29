################FUNCTIONS########

library(tidyverse)

get_pos <- function(x, mafs, fstates, nhaps, N, CC=NULL, type, poslist) {
  
  if(type=="nsnp")
  {
    half.window.size <- (N-1) / 2 # half window size
    #get mafs for window
    Y <- as.matrix(mafs[(x-half.window.size):(x+half.window.size)],ncol=1)
    #get founder data for window
    X <- as.matrix(fstates[(x-half.window.size):(x+half.window.size),],ncol=nhaps)
    
    pos.set <- poslist[(x-half.window.size):(x+half.window.size),]
    
    if(is.null(CC)) 
      {
      Ccov <- NULL
      } else {
        #get coverage for window
        Ccov <- as.matrix(CC[(x-half.window.size):(x+half.window.size)],ncol=1)
        }
  }else{
    if(type=="ppos")
    {
      #get positions
      xp <- as.numeric(strsplit(rownames(fstates)[x], split="_", fixed=TRUE)[[1]][2])
      in_pos <- which(poslist$Ppos >= (xp - (N/2)) & poslist$Ppos <= (xp + (N/2)))
        
      #get mafs for window
      Y <- as.matrix(mafs[in_pos],ncol=1)
      #get founder data for window
      X <- as.matrix(fstates[in_pos,],ncol=nhaps)
      
      pos.set <- poslist[in_pos,]
      
      if(is.null(CC)) 
      {
        Ccov <- NULL
      } else {
        #get coverage for window
        Ccov <- as.matrix(CC[in_pos],ncol=1)
      }
    }else{
      if(type=="gpos")
      {
        #get positions
        xp <- as.numeric(strsplit(rownames(fstates)[x], split="_", fixed=TRUE)[[1]][2])
        xg <- poslist[which(poslist$Ppos == xp),"Gpos"]
        in_pos <- which(poslist$Gpos >= xg - (N/2) & poslist$Gpos <= xg + (N/2))
        
        #get mafs for window
        Y <- as.matrix(mafs[in_pos],ncol=1)
        #get founder data for window
        X <- as.matrix(fstates[in_pos,],ncol=nhaps)
        
        pos.set <- poslist[in_pos,]
        
        if(is.null(CC)) 
        {
          Ccov <- NULL
        } else {
          #get coverage for window
          Ccov <- as.matrix(CC[in_pos],ncol=1)
        }
      }
    }
  }
  
  if(nrow(pos.set) > 1)
  {
    
    vvs <- apply(X, 2, var)
    
    if(all(vvs > 0)){
      fcor <- cor(X)
      fcor <- fcor[upper.tri(fcor)]
    }else{
      fcor <- NA
    }
    
  }else{
    fcor <- NA
  }
  
  
  
  interval_info <- data.frame(
    "ChArm" = poslist$ChArm[x],
    "Ppos" = poslist$Ppos[x],
    "type" = type,
    "nsnp" = nrow(pos.set),
    "Pdist" = pos.set$Ppos[nrow(pos.set)] - pos.set$Ppos[1],
    "Gdist" = pos.set$Gpos[nrow(pos.set)] - pos.set$Gpos[1],
    "cormax" = max(abs(fcor)),
    "nHcor" = length(which(abs(fcor)>0.75))
  )
  
  return(list(Y = Y, X = X, Ccov = Ccov, IntInfo = interval_info))
}


max.mle <- function(hap.freq, pool.cts, ff)
{
  ex.freq <- colSums(t(ff)*hap.freq)
  ex.freq[ex.freq < 0.0001] <- 0.0001
  ex.freq[ex.freq > 0.9999] <- 0.9999
  LL <- -sum(dbinom(pool.cts[,1],pool.cts[,2],ex.freq, log=TRUE))
  return(LL)
}



est.mle <- function(X,Y, CC=NULL, nhaps = 8, p){
 
  est <- numeric(length=nhaps) #set up output
  
  if(is.null(CC)) {
    Ccov <- rep(100, nrow(X))   
    pool.cts <- cbind(round(as.numeric(Y*Ccov)), Ccov)
 
    out <- optim(par=p, max.mle, pool.cts=pool.cts, 
                 ff=X, hessian =TRUE, method = "L-BFGS-B",
                 lower = rep(0, nhaps), upper = rep(1, nhaps))
      
  } else {
    Ccov <- CC
    pool.cts <- cbind(round(as.numeric(Y*Ccov)), Ccov)
    
    out <- optim(par=p, max.mle, pool.cts=pool.cts, 
                 ff=X, hessian =TRUE, method = "L-BFGS-B",
                 lower = rep(0, nhaps), upper = rep(1, nhaps))
  }
  
  # haplotype frequencies
  est <- out$par
  
  return(as.numeric(est))
} #function close


#function 1: maximizer. This gets a quantity analagous to the 
#sum of the squared differences while also trying to account for coverage

max.optim <- function(p,X,Y,N,Ccov=NULL,w) {
  # p is set of estimated haplotype frequencies. [Vector that sums to one, all equal.]
  # X is the set of haplotype allele frequencies (Vector of 8/15 1s or 0s). 
  # Y is the set of allele frequencies from pooled sample
  # N is window size
  # Ccov is coverage
  if(is.null(Ccov)) {
    SS1 <- sum((Y - X %*% p)^2) # SS = Sum of Squares
  } else {
    SS1 <- sum(sqrt(Ccov)*((Y - X %*% p)^2)) # SS = Sum of Squares
  }
  #w <- 100 #serves as a strong penalty if sum of types > 1  
  SS2 <- w * ((sum(p) - 1)^2) 
  return(SS1 + SS2)
}


est.optim <- function(X, Y, CC = NULL, nhaps = 8, pen = 100, p) {
  out <- optim(p, max.optim, X = X,Y = Y, N = N, w = pen, Ccov = CC,
               method = "L-BFGS-B",
               lower = rep(0, nhaps), upper = rep(0.95, nhaps))
  return(as.numeric(out$par))
}



get_hap_freqs <- function(probs, pool)
{
  #probs = haplotype prob data at position with one column labeled RIL. all other columns 
  #should be haplotype probability columns
  
  #pool = vector of RILs in pool
  
  prob.set <- probs[as.character(pool),]
  #prob.set <- prob.set[,-which(colnames(prob.set)=='RIL')]
  return(colSums(prob.set)/nrow(prob.set))
  
}

infer.pos.geno<-function(chrom, newpos, poslist, genolist)
{
  #chrom = chromosome
  #newpos = any position
  #poslist = position list from DSPR library
  #genolist = complete list of haplotype probabilities in same order as poslist
  
  lowpos<-as.integer(newpos/1e4)*1e4
  highpos<-as.integer((newpos+1e4)/1e4)*1e4
  
  if(lowpos < min(poslist[poslist$chr==chrom,'Ppos']))
  {
    lowpos <- min(poslist[poslist$chr==chrom,'Ppos'])
    highpos <- min(poslist[poslist$chr==chrom,'Ppos'])+10000
    
  }
  
  if(highpos > max(poslist[poslist$chr==chrom,'Ppos']))
  {
    lowpos <- max(poslist[poslist$chr==chrom,'Ppos'])-10000
    highpos <- max(poslist[poslist$chr==chrom,'Ppos'])
    
  }
  
  low.dat <- genolist[[which(poslist$chr==chrom & poslist$Ppos==lowpos)]]
  high.dat <- genolist[[which(poslist$chr==chrom & poslist$Ppos==highpos)]]
  
  slope<-(high.dat-low.dat)/(highpos-lowpos)
  intercept<- -slope*highpos+high.dat
  newprobs<-slope*newpos+intercept 
  #newprobs$RIL <- rillistl
  return(newprobs)
  
}

df_to_stan_data <- function(Y, M) {
  list (Y = as.numeric(Y),
        V1 = M[, "freq_A1"],
        V2 = M[, "freq_A2"],
        V3 = M[, "freq_A3"],
        V4 = M[, "freq_A4"],
        V5 = M[, "freq_A5"],
        V6 = M[, "freq_A6"],
        V7 = M[, "freq_A7"],
        V8 = M[, "freq_AB8"],
        N = nrow(M))
}

precision <- function(s) {return(1 / (s ^ 2))}
