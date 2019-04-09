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

robustify.scale <- function(x, w = 9, center = 0) (biweight.scale.est(x = x, center = center, w = w)/sqrt(mean((x-center)^2)))


FAST <- function(spm, method = "robust",mask = NULL, alpha = 0.05, 
                 verbose=FALSE,all=FALSE, stopping=TRUE, K = 100)
{
    ny <- dim(spm) ## size of image
    n <- prod(ny) ## number of voxels
    ##************************************
    ## create mask if non was provided
    ##***********************************
    if(is.null(mask)){
        maskline <- apply(abs(spm), 1:length(dim(spm)), median)
        mask <- ifelse(maskline > quantile(abs(spm), probs = 0.50), TRUE, FALSE)
        dim(mask) <- ny
    }

    if(is.null(mask)){
        mask <- array(TRUE,dim=ny)
    }
    ## choose between am-fast and ar-robust
    method <- match.arg(method,choices = c("model-based","robust","likelihood")) 
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
    lny <- length(ny)
    iotaG <- qgumbel(p=1-alpha) ## k=1, Gumbel(x)
    iotaW <- qrweibull(p=1-alpha)  ## After k>1, limiting distributionis  reverse weibull 
    ny2 <- length(dim(Zmap))
    
    if(verbose){
        cat(paste("",paste(rep("-",100),collapse=""),"\n",sep="")) 
        cat("\t | k | a^{(k)}_n | b^{(k)}_n |   FWHM^{(k)}   | rho^{(k)} | tau^{(k)} |\n")
        cat(paste("",paste(rep("-",100),collapse=""),"\n",sep="")) 
    }
    ##********************************************
    ## If not initial values are not provided for gcv
    ## find optimal using gcv 
    ##**************************************
    
    for(k in 1:K){
        if (method == "robust") {
            if (k == 1) {
                ll.fwh.current <- Inf
                for(fwhm.init in seq(from = 0.05, to = 8.1, by = 2)) {
                    ff <- rep(fwhm.init, ny2)
                    gcv.init.est <- gcv.smooth3d.general(y=Zmap,initval=ff)
                    tmp.val <- gcv.init.est$par.val$value
                    if(tmp.val < ll.fwh.current) {
                        ll.fwh.current <- tmp.val
                        gcv.init <- gcv.init.est$par.val$par
                    }
                }
            } else 
                gcv.init <- gcv$par.val$par
            gcv <- gcv.smooth3d.general(y=Zmap,initval=gcv.init)
            Zim <- gcv$im.smooth
        } else {
            if (method == "model-based") {
                if (k == 1) {
                    ll.fwh.current <- Inf
                    for(fwhm.init in seq(from = 0.05, to = 1, by = 0.15)) {
                        ff <- rep(fwhm.init, ny2+1)
                        ff[1] <- ff[1]*10
                        mb.init.est <- mb.smooth3d.general(y=Zmap,initval=ff)
                        tmp.val <- mb.init.est$par.val$value
                        if(tmp.val < ll.fwh.current) {
                            ll.fwh.current <- tmp.val
                            mb.init <- mb.init.est$par.val$par
                        }
                    }
                } else 
                    mb.init <- c(1,mb$par.val$par[-1])
                mb <- mb.smooth3d.general(y=Zmap,initval=mb.init)
                Zim <- mb$im.smooth
            } else
                Zim <- Zmap
        }
        ##*********************************************
        ## obtain the FWHM corresponding to (7) of paper, 
        ## this will be used in the AM-FAST smoothing and 
        ## also in the AR-FAST analysis
        ##*********************************************
        
        if (k == 1) {
            llh.fwh.current <- -Inf
            llhd.est <- list(value = Inf, par = c(1,rep(0.01, ny2)))
            for(ff in seq(from=1,to=8,by = 2)){
                llhd.est2 <-try(optim(par = c(1,rep(ff,ny2)), fn = fwhm2.llhd.wrapper,
                                      tstat=Zim, eps = 1e-16,
                                      control=list(fnscale=-1)),silent=TRUE)
                if(class(llhd.est2) != "try-error"){
                    tmp.val <- llhd.est2$value
                    if(tmp.val > llh.fwh.current){
                        llh.fwh.current <- tmp.val
                        llhd.est <- llhd.est2
                    }
                }
            }
        }else {
            est.par <- c(1, FWHM[k-1,])
            llhd.est <- optim(par = est.par, fn = fwhm2.llhd.wrapper,
                              tstat=Zim, eps = 1e-16, control=list(fnscale=-1))
        }

#        browser()
        
        ## ************************************************
        ## This bandwidth is the one used in AR-FAST Thresholding.
        ## It is the bandwidth used in the AM-FAST Smoothing 
        ## in lieu of a direct GCV estimate of h
        ##***************************************************

        if(method == "model-based"){
            ##****************************************
            ## We need to smooth current iteration's 
            ## beginning SPM with Gaussian with FWHM h
            ## Compute T^*_{h_k} smoothing Z-map using optimal fhwm
            ##***************************************************
            mb$im.smooth <- mb$im.smooth/(llhd.est$par[1]*robustify.scale(mb$im.smooth[mask]))
            Zmap <- mb$im.smooth ##Smooth Map under robust smoothing
            fwhm.est <- llhd.est$par[-1] ## Estimated FWHM  
            FWHM[k,] <- fwhm.est
            ## LogLikelihood value
            logLike[k] <- fwhm.llhd.wrapper(tstat = Zmap,fwhm=fwhm.est)
        } else {        
            if(method == "robust"){
                gcv$im.smooth <- gcv$im.smooth/(llhd.est$par[1]*robustify.scale(gcv$im.smooth[mask]))
                Zmap <- gcv$im.smooth ##Smooth Map under robust smoothing
                fwhm.est <- llhd.est$par[-1] ## Estimated FWHM  
                FWHM[k,] <- fwhm.est
                ## LogLikelihood value
                logLike[k] <- fwhm.llhd.wrapper(tstat = Zmap,fwhm=fwhm.est)
            } else {
                ##****************************************
                ## We need to smooth current iteration's 
                ## beginning SPM with Gaussian with FWHM h
                ## Compute T^*_{h_k} smoothing Z-map using optimal fhwm
                ##***************************************************
                Zmap <- Gauss.smooth(tstat=Zmap,fwhm=llhd.est$par[-1])
                llhd3.est <- optim(par = llhd.est$par, fn = fwhm2.llhd.wrapper,
                                   tstat=Zmap, eps = 1e-16,
                                   control = list(fnscale=-1))
                Zmap <- Zmap/(llhd3.est$par[1]*robustify.scale(Zmap[mask]))
                fwhm.est <- llhd.est$par[-1] ## Estimated FHWM
                FWHM[k,] <- fwhm.est  
                logLike[k] <- fwhm.llhd.wrapper(tstat = Zmap,fwhm=fwhm.est)
            }
        }
        ## compute varrho_k = R^(1/2) 1
        varrho[k] <- var.rho(n=ny,fwhm = fwhm.est)
        if(lny==2){
            Zmap.sm[,,k] <- Zmap
        } else{
            Zmap.sm[,,,k] <- Zmap
        }
        
        if(k == 1){
            n.not.act <- sum(mask)
            bnk[k] <- qnorm(p = 1-1/n.not.act,sd = varrho[k]) ## bn = F^(1-1/n)/rho
            ank[k] <- 1/(n.not.act* dnorm(x = bnk[k],sd = varrho[k]))## an = 1/(n * rho*f(bn))
            tauk[k] <- (ank[k] * iotaG + bnk[k]) ## Threshold at Gumbel Step
            ## Determining which voxels are activated if > tau_k
            zeta[(Zmap > tauk[k]) & (mask)] <- 1
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
            if(sum(mask) == sum(zeta.old)){
                break;
            }
            n.not.act <- sum((zeta.old==0)&(mask)) ## number of voxels inside the RIO that still not activated 
            pp <- 1 - 1/n.not.act
            if((pp <= 0.0)|(pp >= 1.0)){
                k <- k-1
                break;
            }
            ## Truncation at Reverse Weibull
            etak[k-1] <- tauk[k-1] ## Truncation at Reverse Weibull
            bnk[k] <- etak[k-1]
            ank[k] <- (etak[k-1] - qtruncnorm(p=pp,a = -Inf,b = etak[k-1],sigma=varrho[k]))
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

        if  (sum(Zeta[,k])==0)
            break;
        
        if(k >= 3){
            ##*******************************************
            ## Compute Jaccard Index between activation map
            ## at kth step against step (k-1)th step and jaccard
            ## index kth between (k+1)th step
            ##****************************************
            x1 <- Zeta[,k-2]          
            x2 <- Zeta[,k-1] 
            x3 <- Zeta[,k]      
            JI[k-2] <- jaccard.index(x = x1,y = x2)
            JI[k-1] <- jaccard.index(x = x2,y = x3)
            if(is.na(JI[k-2])){
                JI[k-2] <- 0
            }
            
            if(is.na(JI[k-1])){
                JI[k-1] <- 0
            }
            
            if(stopping){
                if ((JI[k-2] >= JI[k-1])|(JI[k-2] >= 0.9)){
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
            Zmap.sm <- Zmap.sm[,,,1]
        } else{
            ActMap.step <- Zeta
            Zeta <- Zeta[,,,dim(Zeta)[4]-1]
            Zmap.sm <- Zmap.sm[,,,1:k]
        }
    } else{
        if(k == 1){
            ActMap.step <- Zeta
            Zeta <- Zeta[,,1]
            Zmap.sm <- Zmap.sm[,,1]
        } else{
            ActMap.step <- Zeta
            Zeta <- Zeta[,,dim(Zeta)[3]-1]
            Zmap.sm <- Zmap.sm[,,1:k]
        }
    }
    if(verbose){
        cat("Total activated voxels:", sum(Zeta==1),"\n")
        cat(paste("##",paste(rep("=",100),collapse=""),"##\n",sep=""))
    }
    if(all){
        list(FWHM=FWHM[1:k,],BW=bw[1:k], ActMap=Zeta, SPM = spm, SmoothSPMStep=Zmap.sm,
             SmoothSPM=Zmap, JaccardIndex=JI[1:(k-1)],TrackMap = ActMap.step,ThresholdStep=tauk[1:k],
             TruncationStep= etak[1:k],VarRho = varrho[1:k],LogLike=logLike[1:k],
             An = ank[1:k], Bn = bnk[1:k])
    } else{
        list(ActMap=Zeta, SPM = spm,SmoothSPM=Zmap)
    }

}
## wrapper function

FASTfMRI <- function(spm,alpha=0.05,method="AR",two.sided=FALSE,...){

    ## AM: model-based adaptive smoothing
    ## AR: robust adaptive smoothing
    ## AL: likelihood adaptive smoothing
    
    method <- match.arg(method,choices = c("AM","AR","AL"))
    method <- ifelse(method=="AM","model-based",ifelse(method=="AL","likelihood","robust"))
    if(two.sided){
        ff1 <- FAST(spm = spm, method = method, alpha = alpha/2,...)
        ff2 <- FAST(spm = -spm, method = method, alpha = alpha/2,...)
        ## compute union for two-sided test only. See Almodovar and Maitra (2018) for more details
        ff <- list()
        ff$ActMap <- as.numeric(ff1$ActMap | ff2$ActMap)
        dim(ff$ActMap) <- dim(spm)
        ff$SmoothMap <- ff1$SmoothSPM*ff1$ActMap + ff2$SmoothSPM*ff2$ActMap
        ff$SPM <- spm
    } else{
        ff <- FAST(spm = spm, method = method, alpha = alpha,...)
    }
    
    ff
}
