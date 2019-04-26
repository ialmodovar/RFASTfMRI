##*********************************************************************
##
## @file: smoothing_fwhm.R
##
## Perform gaussian smoothing and maximum likelihood for AM-FAST
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
##**********************************************************************


gaussian.filter.1d <- function(n, fwhm = 1, scaled= TRUE)
{
    spar <- fwhm/2.3548
    if (spar < 1e-16) {
        ff <- rep(0,n)
        ff[1] <- 1
    } else{
        ff <- exp(-((0:(n/2))/spar)^2/2)
        ff <- c(ff, rev(ff[2: ((n + 1)/2)]))
    }
    if(scaled)
        ff <- ff/sum(ff)
    ff
}

gaussian.filter.2d <- function(n1, n2 = n1, fwhm1, fwhm2 = fwhm1, scaled = TRUE)
{
    outer(gaussian.filter.1d(n = n1, fwhm = fwhm1, scaled = scaled), gaussian.filter.1d(n = n2, fwhm = fwhm2, scaled = scaled))
}

gaussian.filter.3d <- function(n1, n2 = n1, n3 = n1, fwhm1, fwhm2 = fwhm1, fwhm3 = fwhm2, scaled = TRUE)
{
    f <- array(rep(gaussian.filter.1d(n = n3, fwhm = fwhm3, scaled = scaled), each = n1 * n2), dim = c(n1, n2, n3))
    g <- array(gaussian.filter.2d(n1 = n1, n2 = n2, fwhm1 = fwhm1, fwhm2 = fwhm2, scaled = scaled), dim = c(n1, n2, n3))
    f * g    
}

scaled.fft <- function(z, inverse = FALSE) (fft(z = z, inverse = inverse)/sqrt(ifelse(is.null(dim(z)), length(z), prod(dim(z)))))

ffthalf.gaussian <- function (n, fwhm, scaled = FALSE, eps = 0.01)
{
    if (length(fwhm) == 2)
        f <- gaussian.filter.2d(n1 = n[1], n2 = n[2], fwhm1 = fwhm[1],
            fwhm2 = fwhm[2], scaled = scaled)
    else {
        if (length(fwhm) == 3)
            f <- gaussian.filter.3d(n1 = n[1], n2 = n[2], n3 = n[3],
                fwhm1 = fwhm[1], fwhm2 = fwhm[2], fwhm3 = fwhm[3],
                scaled = scaled)
        else stop("this currently only works for 2d or 3d")
    }
    f.fft <- scaled.fft(f)
    f.1 <- sqrt(f.fft/sqrt(prod(n)))
    Re(f.1)
}


 fftminushalf.gaussian <- function (n, fwhm, scaled = FALSE, eps = 0.01)
{
   f.1 <- ffthalf.gaussian(n = n, fwhm = fwhm, scaled = scaled)
   f.2 <- 1/f.1
   f.2[abs(f.1) < eps * max(abs(f.1))] <- 0
   f.2/sqrt(prod(n))
}

gaussian.minushalf <- function(n, fwhm, scaled = FALSE, eps = 0.01) (Re(scaled.fft(fftminushalf.gaussian(n = n, fwhm = fwhm, scaled = scaled, eps = eps), inverse = TRUE))/sqrt(prod(n)))

gaussian.half <- function(n, fwhm, scaled = FALSE, eps = 0.01) (Re(scaled.fft(ffthalf.gaussian(n = n, fwhm = fwhm, scaled = scaled, eps = eps), inverse = TRUE)))

Gauss.smooth <- function(tstat, fwhm,scaled=TRUE) {
    ##
    ## tstat = 2- or 3-d test statistics
    ## fwhm = vector of fwhms of length 2 or 3.
    ##
    n <- dim(tstat)
    sum.t <- sum(tstat)
    if (length(fwhm) == 2)
        f <- gaussian.filter.2d(n1 = n[1], n2 = n[2], fwhm1 = fwhm[1], fwhm2 = fwhm[2], scaled = scaled)
    else {
        if (length(fwhm) == 3)
            f <- gaussian.filter.3d(n1 = n[1], n2 = n[2], n3 = n[3], fwhm1 = fwhm[1], fwhm2 = fwhm[2], fwhm3 = fwhm[3], scaled = scaled)
        else
            stop("only works for 2d or 3d")
    }
    f.fft <- fft(z = f)
    t.fft <- fft(z = tstat)
    tf.fft <- t.fft * f.fft
    tf <- Re(fft(z = tf.fft, inverse = TRUE))
    sum.tf <- sum(tf)
    tf <- tf * sum.t/sum.tf
}
##***************
## L(FWHM|data)
##**************

fwhm.llhd <- function(fwhm, tstat, eps = 1e-16) {
    ## 
    ## find the decorrelated test statistics using the Gaussian kernel of 
    ## FWHM fwhm
    ##
    
    f1 <- fftminushalf.gaussian(n = dim(tstat), fwhm = fwhm, scaled = FALSE, eps = eps)
    X.cor <- Re(scaled.fft(scaled.fft(tstat)*f1, inverse = TRUE))
    (-sum(X.cor^2)/2 + sum(log((f1[f1 > 0]))))
}

##****************
## Restrict h_ set
##**************
var.rho <-  function(n, fwhm, eps = 1e-16) {
    ##
    ## n = dim of the image
    ## fwhm = fwhms of the final smoother
    ##
      sum(gaussian.half(n, fwhm = fwhm, eps = eps))
}


fwhm.llhd.wrapper <- function(fwhm, tstat, eps = 1e-16) (ifelse( ((min(fwhm) < 1e-10)), -Inf,  fwhm.llhd(fwhm = fwhm, tstat = tstat, eps = eps)))

fwhm2.llhd <- function(pars, tstat, eps = 1e-16) {
    ## 
    ## find the decorrelated test statistics using the Gaussian kernel of 
    ## FWHM fwhm
    ##

    sigma <- pars[1]
    fwhm <- pars[-1] 
    
    f1 <- fftminushalf.gaussian(n = dim(tstat), fwhm = fwhm, scaled = FALSE, eps = eps)
    X.cor <- Re(scaled.fft(scaled.fft(tstat)*f1, inverse = TRUE))
    (-sum(X.cor^2)/(2*sigma^2) + sum(log((f1[f1 > 0]))) -log(sigma)*sum(f1>0))
}

diff.pars <- function(x,...) (outer(x,x,FUN=function(x,y)(abs(x-y)),...) )

fwhm2.llhd.wrapper <- function(pars, tstat, eps = 1e-16) (ifelse( ((pars[1]<=0) | (min(pars[-1]) < 0.5) | (max(pars[-1]) > 8) | (max(diff.pars(pars[-1])) > 2)), -Inf,  fwhm2.llhd(pars = pars, tstat = tstat, eps = eps)))

