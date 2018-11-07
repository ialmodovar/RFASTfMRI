##====================================================================
##
## @file: smooth3d.R
##
## Perform robust smoothing of Garcia 2010 for 3D dataset.
##
## Require: fftw
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
## Author:
##
## Ranjan Maitra, 
## Department of Statistics
## Iowa State University
## maitra@iastate.edu
## Modified by
## Israel A Almodovar-Rivera, 
## Department of Biostatistics and Epidemiology
## Graduate School of Public Health
## University of Puerto Rico at Medical Science Campus
## israel.almodovar@upr.edu 
##====================================================================




DCT2 <- function(x, inverse = FALSE, type=2)
  {
# from
# https://stackoverflow.com/questions/11215162/how-to-perform-a-fast-dct-discrete-cosine-transform-in-r
#
#
    require(fftw)
    if (is.vector(x)) {
      if (inverse)
        IDCT(x, type = type)
      else 
        DCT(x, type = type)
    }
    else {
      y <- x
      if (inverse) {
        y <- t(apply(X = y, MAR = 1, FUN = IDCT, type = type))
        y <- apply(X = y, MAR = 2, FUN = IDCT, type = type)
      }
      else {
        y <- t(apply(X = y, MAR = 1, FUN = DCT, type = type))
        y <- apply(X = y, MAR = 2, FUN = DCT, type = type)
      } 
      y      
    }
  }

DCT3 <- function(x, inverse = FALSE, type = 2)
  {
    if (length(dim(x)) <=2) 
      DCT2(x, inverse = inverse, type = type)
    else {
      if (length(dim(x)) > 3)
        cat("DCT not implemented for arrays of dimension more than 4\n")
      else {
        y <- x
        y <- array(apply(X = y, MAR = 3, FUN = DCT2, inverse = inverse, type = type), dim = dim(x))
        y <- aperm(apply(X = y, MAR = c(1, 2), FUN = DCT2, inverse = inverse, type = type), perm = c(2, 3, 1))
        y
      } 
    }
  }


gcv.score <- function(s, Lambda, DCT3y, N)
{
    if (s < 0)
        Inf
    else {
        Gamma <- 1/(1 + s * Lambda^2)
        RSS <- sum((DCT3y * (Gamma - 1))^2)
        TrH <- sum(Gamma)
        GCVs <- RSS / prod(N) / (1 - TrH/prod(N))^2
                                        #     print(c(s,GCVs))
        GCVs
    }
}


gcv.smooth3d <- function(y, interval)
{
    if (is.null(dim(y)))
        n <- length(y)
    else
        n <- dim(y)
    if (length(n) == 1) 
        lambda <- -2 + 2 * cos((0:(n - 1)) * pi / n)
    else {
        if (length(n) == 2) {
            lambda <- kronecker( X = -2 + 2 * cos((0:(n[1] - 1)) * pi / n[1]),
                                Y = -2 + 2 * t(cos((0:(n[2] - 1)) * pi / n[2])), 
                                FUN = "+")
        }
        if (length(n) == 3) {
            lambda <- array(kronecker( X = -2 + 2 * cos((0:(n[1] - 1)) * pi / n[1]),
                                      Y = -2 + 2 * t(cos((0:(n[2] - 1)) * pi / n[2])), 
                                      FUN = "+"), dim = n) +
                array(rep(-2 + 2 * cos((0:(n[3] - 1)) * pi / n[3]), each = n[1]*n[2]), dim = n)
        }
    }
    dct3y <- DCT3(y, inverse = FALSE)
    
    par.val <- optimize(gcv.score, interval = interval, Lambda = lambda,DCT3y = dct3y, N = n, tol = .Machine$double.eps)
    
    shat <- par.val$minimum

    gamma <- 1/(1 + shat * lambda^2)

    c(list(im.smooth = DCT3(gamma * dct3y, inverse = TRUE), par.val = par.val))
}


setup.eigvals <- function(s, n)
{
    if (length(n) == 2) {
        Lambda <- kronecker( X = (-2 + 2 * cos((0:(n[1] - 1)) * pi / n[1]))*s[1],
                            Y = (-2 + 2 * t(cos((0:(n[2] - 1)) * pi / n[2])))*s[2], 
                            FUN = "+")
    }
    else {
        if (length(n) == 3) {
            Lambda <- array(kronecker( X = (-2 + 2 * cos((0:(n[1] - 1)) * pi / n[1]))*s[1],
                                      Y = (-2 + 2 * t(cos((0:(n[2] - 1)) * pi / n[2])))*s[2], 
                                      FUN = "+"), dim = n) +
                array(s[3]*rep(-2 + 2 * cos((0:(n[3] - 1)) * pi / n[3]), each = n[1]*n[2]), dim = n)
        }
    }
    Lambda
}

gcv.score.general <- function(s, DCT3y)
{
    if (min(s) < 0)
        Inf
    else {
        n <- dim(DCT3y)
        Lambda <- setup.eigvals(s, n)
        Gamma <- 1/(1 + Lambda^2)
        RSS <- sum((DCT3y * (Gamma - 1))^2)
        TrH <- sum(Gamma)
        GCVs <- RSS / prod(n) / (1 - TrH/prod(n))^2
        GCVs
    }
}


gcv.smooth3d.general <- function(y, initval)
{
    if (length(dim(y)) != length(initval))
        cat("Error: use of this function only makes sense for higher dimensions, use gcv.smooth3d instead\n")
    else {
        n <- dim(y)
        dct3y <- DCT3(y, inverse = FALSE)
        
        par.val <- optim(par = initval, fn = gcv.score.general, DCT3y = dct3y)
        
        shat <- par.val$par

        lambda <- setup.eigvals(shat, n)
        gamma <- 1/(1 + lambda^2)
        c(list(im.smooth = DCT3(gamma * dct3y, inverse = TRUE), par.val = par.val))
    }
}








setup.eigvals.spatial.only <- function(s, n, MAR)
{
                                        # MAR denotes the index which will not be smoothed (the time series)
                                        #

    if (MAR == 3) {
        Lambda <- array(kronecker( X = (-2 + 2 * cos((0:(n[1] - 1)) * pi / n[1]))*s[1],
                                  Y = (-2 + 2 * t(cos((0:(n[2] - 1)) * pi / n[2])))*s[2], 
                                  FUN = "+"), dim = n)
    }
    else {
        if (MAR == 1) {
            Lambda <- array(kronecker( X = rep(0, n[1]),
                                      Y = (-2 + 2 * t(cos((0:(n[2] - 1)) * pi / n[2])))*s[1], 
                                      FUN = "+"), dim = n) +
                array(s[2]*rep(-2 + 2 * cos((0:(n[3] - 1)) * pi / n[3]), each = n[1]*n[2]), dim = n)
        }
    }
    Lambda
}


gcv.score.general.spatial.only <- function(s, DCT3y, MAR)
{
    if (min(s) < 0)
        Inf
    else {
        n <- dim(DCT3y)
        Lambda <- setup.eigvals.spatial.only(s, n, MAR)
        Gamma <- 1/(1 + Lambda^2)
        RSS <- sum((DCT3y * (Gamma - 1))^2)
        TrH <- sum(Gamma)
        GCVs <- RSS / prod(n) / (1 - TrH/prod(n))^2
        GCVs
    }
}


gcv.smooth3d.general.spatial.only <- function(y, initval, MAR)
{
    n <- dim(y)
    dct3y <- DCT3(y, inverse = FALSE)
    
    par.val <- optim(par = initval, fn = gcv.score.general.spatial.only, DCT3y = dct3y, MAR = MAR)
    
    shat <- par.val$par
    
    lambda <- setup.eigvals.spatial.only(shat, n, MAR)
    gamma <- 1/(1 + lambda^2)
    c(list(im.smooth = DCT3(gamma * dct3y, inverse = TRUE), par.val = par.val))
}


