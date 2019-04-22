\name{FASTfMRI}
\alias{FASTfMRI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Adaptive Smoothing and Thresholding (Wrapper)
}
\description{
Compute FAST for activation detection in fMRI studies. 
}
\usage{
FASTfMRI(spm, alpha = 0.05, method = "AR", two.sided = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spm}{ Statistical Parametric Map. Array of a 2D or 3D array.}
  \item{method}{If method = "AR" will perform robust smoothing of Garcia
    (2010). If method = "AL" will perform maximum likelihood
    estimate of the optimal FWHM to perform smoothing. If method = "AM" will perform model-based smoothing.}
  \item{alpha}{Significance level, default 0.05}
  \item{two.sided}{Perform a two-sided hypothesis. Default is FALSE. If
    true will also perform FAST in a -SPM with a significance level of
    alpha/2. H0: c'b= 0 versus H1: c'b!= 0}
  \item{\dots}{Further internal arguments for the smoothing algorithm usually. Accept mask and verbose.}
}
\details{
Compute FAST for activation detection in fMRI studies. By default perform one sided test.

FWHM: a 2D or 3D array containing the value of the estimated FWHM at each iteration.

ActMap: a 2D or 3D array of 0 and 1. If a voxel is determined to be activated it is classified as 1.

SPM: initial statistical parametric map.

SmoothSPM: Final smoothed map.

JaccardIndex: Jaccard index similirity coefficient at each iteration.

VarRho: Return R^(1/2) 1 at each step.

An: Non-negative Normalizing  constant of the extreme value theory at
k=1, belong to Gumbel domain, for k> 1 belong to Reverse Weibull.

Bn: Normalizing  constant of the extreme value theory at k=1,
belong to Gumbel domain, for k> 1 belong to Reverse Weibull

}
\references{
Almodóvar-Rivera, I., & Maitra, R. (2017). FAST adaptive smoothing and thresholding for improved activation detection in low-signal fMRI. arXiv preprint arXiv:1702.00111.

Garcia, D. (2010). Robust smoothing of gridded data in one and higher dimensions with missing values. Computational statistics & data analysis, 54(4), 1167-1178.
}
\author{
Israel Almodóvar-Rivera and Ranjan Maitra.
}


\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
 
\examples{
## Blue to Red palette from RColorBrewer
burd <- colorRampPalette(rev(c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7",
"#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3","#2166AC","#053061")))

data(Hoff)
data(tmap)
Kmax <- 50
stp <- TRUE

## Create Hoff Mask

mask <- Hoff
mask[is.na(mask)] <- 0
mask <- as.logical(mask)
dim(mask) <- c(128,128)

Hoff2 <- -!Hoff
Hoff2[is.na(Hoff2)] <- -2
mask <- !is.na(Hoff)
in.mask.rows <- (1:nrow(Hoff2))[X = apply(Hoff2 ==0, MARGIN = 1, FUN = any)]
in.mask.cols <- (1:ncol(Hoff2))[X = apply(Hoff2 ==0, MARGIN = 2, FUN = any)]

## Perform AM-FAST (Model-based)

ff.am <- apply(tmap,3,FASTfMRI,mask = mask,method="AM",alpha=0.05,K = Kmax,all=TRUE,stopping=stp,verbose=TRUE)
par(mfrow=c(3,5),mar=rep(0,4))
for(i in 1:15){
  aa <- ff.am[[i]]$ActMap
  tt <- tmap[,,i]
  tmax <- max(tt,na.rm=TRUE)
    tt[aa==0] <- NA
  image(Hoff[in.mask.rows,in.mask.cols],col=grey(1-c(0,0.05,0.025)),axes=FALSE)
  image(tt[in.mask.rows,in.mask.cols],col=burd(16),zlim=c(-tmax,tmax),axes=FALSE,add=TRUE)
  contour(Hoff2[in.mask.rows, in.mask.cols], drawlabels=FALSE, method = "simple", nlevels = 1, add = TRUE, lwd = 0.5)
}

## Perform AR-FAST (robust)

ff.ar <- apply(tmap,3,FASTfMRI,mask = mask,method="AR",alpha=0.05,K = Kmax,all=TRUE,stopping=stp,verbose=TRUE)

for(i in 1:15){
  aa <- ff.ar[[i]]$ActMap
  tt <- tmap[,,i]
  tmax <- max(tt,na.rm=TRUE)  
  tt[aa==0] <- NA
  image(Hoff[in.mask.rows,in.mask.cols],col=grey(1-c(0,0.05,0.025)),axes=FALSE)
  image(tt[in.mask.rows,in.mask.cols],col=burd(16),zlim=c(-tmax,tmax),axes=FALSE,add=TRUE)
  contour(Hoff2[in.mask.rows, in.mask.cols], drawlabels=FALSE, method = "simple", nlevels = 1, add = TRUE, lwd = 0.5)
}


## Perform AL-FAST (likelihood)

ff.al <- apply(tmap,3,FASTfMRI,mask = mask,method="AL",alpha=0.05,K = Kmax,all=TRUE,stopping=stp,verbose=TRUE)
for(i in 1:15){
  aa <- ff.al[[i]]$ActMap
  tt <- tmap[,,i]
  tmax <- max(tt,na.rm=TRUE)
    tt[aa==0] <- NA
  image(Hoff[in.mask.rows,in.mask.cols],col=grey(1-c(0,0.05,0.025)),axes=FALSE)
  image(tt[in.mask.rows,in.mask.cols],col=burd(16),zlim=c(-tmax,tmax),axes=FALSE,add=TRUE)
  contour(Hoff2[in.mask.rows, in.mask.cols], drawlabels=FALSE, method = "simple", nlevels = 1, add = TRUE, lwd = 0.5)
}

## Display
data(HoffActiv)
HoffActiv[!is.na(HoffActiv)] <- 1
HoffActiv[is.na(HoffActiv)] <- 0

jac.am <- sapply(ff.am,function(z) jaccard.index(HoffActiv,z$ActMap))
jac.ar <- sapply(ff.ar,function(z) jaccard.index(HoffActiv,z$ActMap))
jac.al <- sapply(ff.al,function(z) jaccard.index(HoffActiv,z$ActMap))

JacIndex <- data.frame(JaccardIndex = c(jac.am,jac.ar,jac.al), Method =
c(rep("AM-F",length(jac.am)),rep("AR-F",length(jac.ar)),rep("AL-F",length(jac.al))))


with(JacIndex, boxplot(JaccardIndex~Method,ylim=c(0,1)))

## compute total of activated voxels

act.am <- sapply(ff.am,function(z) sum(z$ActMap))
act.ar <- sapply(ff.ar,function(z) sum(z$ActMap))
act.al <- sapply(ff.al,function(z) sum(z$ActMap))

ActVoxels <- data.frame(ActVx = c(act.am,act.ar,act.al), Method =
c(rep("AM-F",length(act.am)),rep("AR-F",length(act.ar)),rep("AL-F",length(act.al))))

with(ActVoxels, boxplot(ActVx~Method))
abline(h = 138,lwd=2)

}

