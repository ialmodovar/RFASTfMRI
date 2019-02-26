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
  tauk <- rep(0,K) ## normalizing sequences value a_n > 0
  etak <- rep(NA,K) ## normalizing sequences value a_n > 0
  Loglike.iter <- rep(0,K)
  lny <- length(ny)
  if(verbose){
    cat(paste("",paste(rep("-",100),collapse=""),"\n",sep="")) 
    cat("\t | k | a^{(k)}_n | b^{(k)}_n |   FWHM^{(k)}   | rho^{(k)} | tau^{(k)} |\n")
    cat(paste("",paste(rep("-",100),collapse=""),"\n",sep="")) 
  }  
  for(k in 1:K){
      if(k==1){   
          if(is.null(fwhm.init)){
              ll.fwh.current <- -Inf
              ny2 <- length(dim(spm))
              like.init <- rep(0,8)
              cnt <- 1
              for( ff in c(0.1,0.25,0.5,0.75,1,2,3,4)){
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
          varrho[k] <- rho(n=ny,fwhm = fwhm.est)
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
##            ## Obtain Sum of first Row Compute rho = sum r_ij
          varrho[k] <- rho(n=ny,fwhm = fwhm.est)
      }
    
      if(k == 1){
          n.not.act <- sum(mask)
          bnk[k] <- qnorm(p = 1-1/n.not.act) # bn = F^(1-1/n)
          ank[k] <- 1/(varrho[k]* n.not.act * dnorm(x = bnk[k])) # an = 1/(n * rho*f(bn*rho))
          
          ##ank[k] <- (bnk[k]/(1+bnk[k]^2))*varrho[k]
          iotaG <- qgumbel(p=1-alpha) ## gumbel assume P((X_{(n)} - bn)/an < x) -> Gumbel(x)
          tauk[k] <- (ank[k] * iotaG + bnk[k]) ## Threshold at Gumbel Step
          ## Determining which voxels are activated if > Eta_k
          zeta[Zmap > tauk[k]] <- 1
          Zeta[,k] <- zeta
          if(verbose){
              if(lny==3){
                  cat(sprintf("\t | %d | %.6f | %.6f | (%.6f, %.6f, %.6f) | %.6f | %.6f | \n",k,ank[k],bnk[k],FWHM[k,1],FWHM[k,2],FWHM[k,3],varrho[k],tauk[k] ))
              } else{
                  cat(sprintf("\t | %d | %.6f | %.6f | (%.6f, %.6f) | %.6f | %.6f | \n",k,ank[k],bnk[k],FWHM[k,1],FWHM[k,2],varrho[k],tauk[k] ))
              }
          }
      } else {
      zeta.old <- Zeta[,k-1]
      zeta <- zeta.old
      n.not.act <- sum((zeta.old==0)&(mask)) ## number of voxels inside the RIO that still not activated
      if(sum(mask) == sum(zeta.old)){
          break;
      }
      pp <- 1 - 1/n.not.act
      if((pp <= 0.0)|(pp >= 1.0)){
          k <- k-1
        break;
      }
      ## Step 2: use Corollary 2 with the h estimated in Algorithm Step 2 a i and ii to get the
      ## smoothed truncation point. This is \eta_k^\bullet and is the \eta to be used in the
      ## calculation of a_n and b_n in Corollary 5 and thus \eta_k.
      
      ## From Step 3 onwards, we calculate \eta_k^\bullet a bit differently.
      ## We find \alpha_{k-1} the upper tail probability of \eta_{k-1} with regard to the standard normal distribution.
      ## Use this and the h obtained from Algorithm Step 2 a i and ii in Corollary 2 to get a new truncation point (note: use \alpha_{k-1}
      ## here in place of alpha). Call the new truncation point \eta_k^\bullet and is the
      ## \eta to be used in the calculation of a_n and b_n in Corollary 5 and thus \eta_k.

##This gives us a set of objective \eta_k^\bullet.

  ##    if(k <3){
          etak[k-1] <- tauk[k-1]/varrho[k]
   ##   } else{
     ##     etak[k-1] <- tauk[k-1]
   ##    alpha <- qrweibull(pnorm(etak[k-1],lower.tail=FALSE))
   ##   }
            ## Truncation constant sufficient statistics of eta is the maximum
      ##etak[k-1] <- max(Zmap[(mask) & (zeta.old==0)],na.rm=TRUE)
      ##etak[k-1] <- tauk[k-1]
      ##bnk[k] <- tauk[k-1]/varrho[k-1]
      ## choose normalizing constants for Reverse Weibull
      bnk[k] <- etak[k-1]/varrho[k]
      ank[k] <- (etak[k-1]/varrho[k] - qtruncnorm(p=pp,a = -Inf,b = etak[k-1]/varrho[k]))/varrho[k]
        
      ##--- we can also use a truncated  |x| < eta
      iotaW <- qrweibull(p=1-alpha)  ## After k>1, limiting distributionis  reverse weibull 
      tauk[k] <- (ank[k] * iotaW +bnk[k]) ## Threshold at Reverse Weibull
      zeta[(Zmap > tauk[k]) & (zeta.old == 0) & (mask)] <- 1  ## activated voxel 
      Zeta[,k] <- zeta
      if(verbose){
          if(lny==3){
              cat(sprintf("\t | %d | %.6f | %.6f | (%.6f, %.6f, %.6f) | %.6f | %.6f | \n",k,ank[k],bnk[k],FWHM[k,1],FWHM[k,2],FWHM[k,3],varrho[k],tauk[k] ))
          } else{
              cat(sprintf("\t | %d | %.6f | %.6f | (%.6f, %.6f) | %.6f | %.6f | \n",k,ank[k],bnk[k],FWHM[k,1],FWHM[k,2],varrho[k],tauk[k] ))
          }
      }
      
    }

    ##if(k == 2){
##  m1 <- sum(Zeta[,k-1])
    ##  m2 <- sum(Zeta[,k])
    ##if((m1==0) & (m2==0)){
    ##    break;
    ##  }

 ##   }
    #----------------------
      if(k >= 3){
          ## Compute Jaccard Indexes between activation map
          ## at kth step against step (k-1)th step and jaccard
          ## index kth between (k+1)th step
      x1 <- Zeta[,k-2]
      x2 <- Zeta[,k-1] 
      x3 <- Zeta[,k]      
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
    cat("Total activated voxels:", sum(Zeta==1),"\n")
    cat(paste("##",paste(rep("=",100),collapse=""),"##\n",sep=""))
  }
  if(all){
  list(FWHM=FWHM[1:k,],BW=bw[1:k], ActMap=Zeta, SPM = spm, SmoothSPMStep=Zmap.sm,
       SmoothSPM=Zmap, JaccardIndex=JI[1:(k-1)],TrackMap = ActMap.step,ThresholdStep=tauk[1:k],
       TruncationStep= etak[1:k],Rho = varrho[1:k],LogLike=logLike[1:k],LogLikeIter = Loglike.iter,
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
