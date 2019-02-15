##***********************************************************************
##
## @file: fastfMRI.R
##
## New method to obtain detect activation in fMRI studies
## while performing adaptive smoothing and 
## taking in account the spatial correlation.
## The SPM is assumed to arrive from a Normal field.
##
## Note the first iteration, F(X) belong to Gumbel domain,
## so a_n and b_n are chosen such as   
## n are not activated, initially all voxels are not activated.
##
## After the first iteration F(X) belong to Reverse Weibull domain.
## This is due since we are assuming X_1,\ldots, X_n \sim N(0,1) 
## truncated at some \eta > 0. This lead to a Weibull x <= \eta.
## so a_n and b_n are chosen such as   
## a_n = F^{-1}_X(1) - F^{-1}_X(1-1/n) =  \eta-F^{-1}_X(1-1/n) ; 
## b_n = \eta, n are not activated, initially all voxels are not activated.
##
## If Jaccard Index (JI) JI_k = JI_{k-1} or JI_k < JI_{k-1}) at 
## iteration k haven't change  stop and use zeta map
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
##
## Authors:
## Israel A Almodovar-Rivera, 
## Department of Biostatistics and Epidemiology
## Graduate School of Public Health
## University of Puerto Rico
## Medical Science Campus
## israel.almodovar@upr.edu 
##
## Ranjan Maitra, 
## Department of Statistics
## Iowa State University
## maitra@iastate.edu 
##***********************************************************************

##****************************************
## truncated normal distribution quantile
##***************************************
qtruncnorm <- function(p, a = -Inf, b = Inf, mu = 0,sigma = 1){
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 1) 
    stop(" 'p' must contain probabilities in (0,1)")
  
  alpha <- (a - mu)/sigma
  beta <- (b - mu)/sigma
  sigma * qnorm(p*(pnorm(beta) - pnorm(alpha))+pnorm(alpha)) + mu
}
##****************************************
## Gumbel distribution quantile
##***************************************

qgumbel <- function(p, mu = 0,beta = 1){
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 1) 
    stop(" 'p' must contain probabilities in (0,1)")
  
  mu - beta * log(-log(p))
}
##****************************************
## Reversed Weibull distribution quantile
##***************************************

qrweibull <- function (p, mu = 0, beta = 1, tau = 1, lower.tail = TRUE) 
{
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 1) 
    stop("'p' must contain probabilities in (0,1)")
  
  if (min(beta) < 0 || min(tau) <= 0) 
    stop("invalid arguments")
  
  if (!lower.tail) 
    p <- 1 - p
  
  mu - beta * (-log(p))^(1/tau)
}
##************************************************
## compute jaccard index  similarity coefficient
##************************************************
jaccard.index <- function(x, y) {
  if ((length(x) == length(y)) & (sum(!is.na(x)) == sum(!is.na(y)))) {
    num <- sum((x == 1) & (y == 1), na.rm = T)
    num/((sum(x == 1, na.rm = T) + sum(y == 1, na.rm = T) - num))
  }
  else {
    stop("Error: number of NAs in x != number of NAs in y")
  }
}

FAST <- function(spm, method = "robust",mask = NULL, alpha = 0.05, fwhm.init=NULL,
                 verbose=FALSE,all=FALSE, stopping=TRUE, K = 100)
{
  ny <- dim(spm) ## size of image
  n <- prod(ny) ## number of voxels
  ##************************************
  ## create mask if non was provided
  ##***********************************
  if(is.null(mask) & (length(dim(spm))==3) ){
    maskline <- apply(spm, 1:3, mean)
    mask <- ifelse(maskline > quantile(spm, probs = 0.80), TRUE, FALSE)
    dim(mask) <- ny
  }

  if(is.null(mask)){
    mask <- array(TRUE,dim=ny)
  }
  ## choose between am-fast and ar-robust
  method <- match.arg(method,choices = c("likelihood","robust")) 
  if(verbose){
    cat(paste("##",paste(rep("=",100),collapse=""),"##\n",sep="")) # Just to create long "="
    cat(" \t \t \t Fast Adaptive Smoothing and Thresholding (FAST)  \n")  
    cat("Smoothing method:", method, "\n")
  }
  logLike <- rep(0,K) ## store likelihood values
  bw <- rep(0,K) ## bandwidth for robust
  Zeta <- array(0,dim=c(n,K)) ## activation map
  zeta <- rep(0,n) ## current vector for activation map
  JI <- rep(0,K) ## jaccard index
  varrho <- rep(1,K) ## image correlation
  ank <- rep(0,K) ## normalizing sequences value a_n > 0
  bnk <- rep(0,K) ## normalizing sequences value b_n \in R
  FWHM <- array(0,dim=c(K,length(ny))) ## optimal FWHM 
  Zmap <- spm ## Copy original Map on Zmap which it will be the smooth one 
  Zmap.sm <- array(0,dim=c(ny,K))
  etak <- rep(0,K) ## normalizing sequences value a_n > 0
  Loglike.iter <- rep(0,K)
  lny <- length(ny)
  if(verbose){
    cat(paste("",paste(rep("-",100),collapse=""),"\n",sep="")) 
    cat("\t | k | a^{(k)}_n | b^{(k)}_n |   FWHM^{(k)}   | rho^{(k)} | eta^{(k)} |\n")
    cat(paste("",paste(rep("-",100),collapse=""),"\n",sep="")) 
  }  
  for(k in 1:K){
      
      if(k == 1){
          if(is.null(fwhm.init)){
              ll.fwh.current <- -Inf
              ny2 <- length(dim(spm))
              like.init <- rep(0,7)
              cnt <- 1
              for( ff in c(0.1,0.25,0.5,0.75,1,2,3)){
                  fwhm.init2 <- rep(ff,ny2)
            ll.fwh <- fwhm.llhd.wrapper(fwhm = fwhm.init2,tstat = Zmap)
                  like.init[cnt] <- ll.fwh
                  cnt <- cnt+1
                  if(ll.fwh >= ll.fwh.current){
                      ll.fwh.current <- ll.fwh
                      fwhm.init <- fwhm.init2
                  }
              }
          }
      }
      if(method == "likelihood"){
        llhd.est <- optim(par = fwhm.init, fn = fwhm.llhd.wrapper,
                          tstat=Zmap, eps = 1e-16, control=list(fnscale=-1))

      fwhm.est <- llhd.est$par ## Estimated FHWM
      logLike[k] <- llhd.est$value ## LogLikelihood value
      FWHM[k,] <- fwhm.est  
      Loglike.iter[k] <- logLike[k]
      ## Compute T^*_{h_k} smoothing Z-map using optimal fhwm
      Zmap <- Gauss.smooth(tstat=Zmap, fwhm = fwhm.est)
      if(lny==2){
      Zmap.sm[,,k] <- Zmap
      } else{
        Zmap.sm[,,,k] <- Zmap
      }
        ## Obtain Sum of first Row Compute rho = sum r_ij
        varrho[k] <- rho.am(n=ny,fwhm = fwhm.est)
      }
      
      if(method == "robust"){
          zmap.gcv <- gcv.smooth3d.general(y=Zmap,initval=fwhm.init)
          Zmap <- zmap.gcv$im.smooth ## Smooth map using robust method
          fwhm.est <- zmap.gcv$par.val$par ##estimate fwhm est
          FWHM[k,] <- fwhm.est
          logLike[k] <- fwhm.llhd.wrapper(fwhm = fwhm.est,tstat = Zmap)
          Loglike.iter[k] <- logLike[k]
      if(lny==2){
        Zmap.sm[,,k] <- Zmap
      } else{
        Zmap.sm[,,,k] <- Zmap
      }
            ## Obtain Sum of first Row Compute rho = sum r_ij
        varrho[k] <- rho.ar(n=ny,fwhm = fwhm.est)
    }
    
    if(k == 1){
      n.not.act <- sum(mask)
      bn <- qnorm(p = 1-1/n.not.act)/varrho[k] # bn = F^(1-1/n)/rho
      an <- 1/(varrho[k]* n.not.act * dnorm(x = bn*varrho[k])) # an = 1/(n * rho*f(bn*rho)) 
      ank[k] <- an
      bnk[k] <- bn
      if(verbose){
          if(lny==3){
              cat(sprintf("\t | %d | %.6f | %.6f | (%.6f, %.6f, %.6f) | %.6f | %.6f | \n",k,an,bn,FWHM[k,1],FWHM[k,2],FWHM[k,3],varrho[k],max(Zmap) ))
          } else{
              cat(sprintf("\t | %d | %.6f | %.6f | (%.6f, %.6f) | %.6f | %.6f | \n",k,an,bn,FWHM[k,1],FWHM[k,2],varrho[k],max(Zmap) ))
          }
      }
      zeta <- sign(Zmap * mask)
      iota <- qgumbel(p=1-alpha,mu=0) ## gumbel assume P((X_{(n)} - bn)/an < x) -> Gumbel(x)
      etak[k] <- ank[k] * iota + bnk[k] ## Threshold at Gumbel Step
      ## Determining which voxels are activated if > Eta_k
      zeta[Zmap < etak[k]] <- 0
      Zeta[,k] <- zeta
    } else {
      zeta.old <- Zeta[,k-1]
      zeta <- zeta.old
      n.not.act <- sum((zeta.old==0)&(mask)) ## number of voxels inside the RIO that still not activated
      eta <- max(Zmap[(zeta.old==0)&(mask)]) ## maximum of non-activated voxels that are in ROI
##      eta <- max((Zmap[(zeta.old==0)&(mask)]),1,na.rm=TRUE) ## maximum of non-activated voxels that are in ROI
      pp <- 1 - 1/n.not.act
      if((eta == Inf)|(eta == -Inf)|(pp <= 0.0)|(pp >= 1.0)){
        break;
      }
      ## choose normalizing constants for Reverse Weibull

      bn <- eta/(varrho[k])^2
      an <- bn - qtruncnorm(p=pp,a = -Inf,b = eta/varrho[k])/varrho[k]

      ank[k] <- an
      bnk[k] <- bn
  
      if(verbose){
          if(lny==3){
              cat(sprintf("\t | %d | %.6f | %.6f | (%.6f, %.6f, %.6f) | %.6f | %.6f | \n",k,an,bn,FWHM[k,1],FWHM[k,2],FWHM[k,3],varrho[k],max(Zmap) ))
          } else{
              cat(sprintf("\t | %d | %.6f | %.6f | (%.6f, %.6f) | %.6f | %.6f | \n",k,an,bn,FWHM[k,1],FWHM[k,2],varrho[k],max(Zmap) ))
          }
      }
      
      ##--- we can also use a truncated  |x| < eta
      iota <- qrweibull(p=1-alpha,mu=eta)  ## After k>1, limiting distributionis  reverse weibull 
      etak[k] <- ank[k] * iota +bnk[k] ## Threshold at Reverse Weibull
      zeta[(Zmap > etak[k]) & (zeta.old == 0) & (mask)] <- 1
      Zeta[,k] <- zeta
      }
    #----------------------
    if(k >= 3){
      x1 <- Zeta[,k-2]
      x2 <- Zeta[,k-1] 
      x3 <- Zeta[,k]
      x1[x1==-1] <- 1 # note that we can have negative activations on zeta (-1), JI take only 0 and 1
      # so we replace the negative activated with 1 just for JI computation purposes
      x2[x2==-1] <- 1
      x3[x3==-1] <- 1
      JI[k-2] <- jaccard.index(x = x1,y = x2)
      JI[k-1] <- jaccard.index(x = x2,y = x3)
      if(is.na(JI[k-2])){
        JI[k-2] <- 1
      }
      
      if(is.na(JI[k-1])){
        JI[k-1] <- 1
      }
      if(stopping){
      if((JI[k-2] >= JI[k-1])|(k == K)){
        break;
      }
      }
    }
  }
  if(verbose){
    cat(paste("",paste(rep("-",100),collapse=""),"\n",sep="")) 
  }
  zeta <- Zeta[,k-1] ## activated final map
  Zeta <- Zeta[,1:k]
  dim(Zeta) <- c(ny,k)

  if(lny == 3){
      if(k == 1){
          ActMap.step <- Zeta
          Zeta <- Zeta[,,,1]
  } else{
    ActMap.step <- Zeta
    Zeta <- Zeta[,,,dim(Zeta)[4]-1]
  }
  } else{
  if(k == 1){
    ActMap.step <- Zeta
    Zeta <- Zeta[,,1]
  } else{
    ActMap.step <- Zeta
    Zeta <- Zeta[,,dim(Zeta)[3]-1]
  }
  }
  if(verbose){
    cat("Activated voxels:", sum(Zeta==1),"\n")
    cat(paste("##",paste(rep("=",100),collapse=""),"##\n",sep=""))
  }
  if(all){
  list(FWHM=FWHM[1:k,],BW=bw[1:k], ActMap=Zeta, SPM = spm, SmoothSPMStep=Zmap.sm,
       SmoothSPM=Zmap, JaccardIndex=JI[1:(k-1)],TrackMap = ActMap.step,EtaStep=etak,
       Rho = varrho[1:k],LogLike=logLike[1:k],LogLikeIter = Loglike.iter,
       An = ank[1:k], Bn = bnk[1:k])
} else{
    list(ActMap=Zeta, SPM = spm,SmoothSPM=Zmap)
}

  }
## wrapper function

FASTfMRI <- function(spm,alpha=0.05,method="AR",two.sided=FALSE,fwhm.init=NULL,...){
    
    method <- match.arg(method,choices = c("AM","AR"))
    if(method=="AR"){
     if(two.sided){
      ff1 <- FAST(spm = spm, method = "robust", alpha = alpha/2,fwhm.init = fwhm.init,...)
      ff2 <- FAST(spm = -spm, method = "robust", alpha = alpha/2,fwhm.init = fwhm.init,...)
      ## compute union for two-sided test only. See Almodovar and Maitra (2018) for more details
      ff <- list()
      ff$ActMap <- as.numeric(ff1$ActMap | ff2$ActMap)
      dim(ff$ActMap) <- dim(spm)
      ff$SmoothMap <- ff1$SmoothSPM*ff1$ActMap + ff2$SmoothSPM*ff2$ActMap
      ff$SPM <- spm
     } else{
      ff <- FAST(spm = spm, method = "robust", alpha = alpha,...)
     }
    }
   if(method=="AM"){
     if(two.sided){
      ff1 <- FAST(spm = spm, method = "likelihood", alpha = alpha/2,fwhm.init = fwhm.init,...)
      ff2 <- FAST(spm = -spm, method = "likelihood", alpha = alpha/2,fwhm.init = fwhm.init,...)
      ## compute union for two-sided test only. See Almodovar and Maitra (2018) for more details
      ff <- list()
      ff$ActMap <- as.numeric(ff1$ActMap | ff2$ActMap)
      dim(ff$ActMap) <- dim(spm)
      ff$SmoothMap <- ff1$SmoothSPM*ff1$ActMap + ff2$SmoothSPM*ff2$ActMap
      ff$SPM <- spm
      } else{
      ff <- FAST(spm = spm, method = "likelihood", alpha = alpha,fwhm.init = fwhm.init,...)
     }
   }
    ff
}
