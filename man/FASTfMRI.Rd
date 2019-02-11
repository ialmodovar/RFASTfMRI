\name{Fast Adaptive Smoothing and Thresholding (Wrapper)}
\alias{FASTfMRI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Adaptive Smoothing and Thresholding
}
\description{
Compute FAST for activation detection in fMRI studies. 
}
\usage{
FASTfMRI(spm, alpha = 0.05, method = "AR", two.sided = FALSE,fwhm.init=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spm}{ Statistical Parametric Map. Array of a 2D or 3D array.}
  \item{method}{If method = "AR" will perform robust smoothing of Garcia (2010). If method = "AM" will perform maximum estimate of the optimal FWHM to perform smoothing.
}
  \item{alpha}{Significance level, default 0.05}
  \item{verbose}{ If TRUE will print every calculation step. Default is FALSE.}
  \item{two.sided}{Perform a two-sided hypothesis. Default is FALSE. If
    true will also perform FAST in a -SPM with a significance level of
    alpha/2. H0: c'b= 0 versus H1: c'b!= 0}
  \item{fwhm.init}{Default is NULL. If not initial FWHM values has been
    given it will chose based one from a given range and uses it the one
    that maximize the log-likelihood.}
  \item{\dots}{Further internal arguments for the smoothing algorithm usually. Accept mask and verbose.}
}
\details{
Compute FAST for activation detection in fMRI studies. By default perform one sided test.

FWHM: a 2D or 3D array containing the value of the estimated FWHM at each iteration.

ActMap: a 2D or 3D array of 0 and 1. If a voxel is determined to be activated it is classified as 1.

SPM: initial statistical parametric map.

SmoothSPM: Final smoothed map.

JaccardIndex: Jaccard index similirity coefficient at each iteration.

Rho: Return R^(1/2) 1 at each step.

An: Normalizing non-negative constant of the extreme value theory.

Bn: Normalizing constant of the extreme value theory.

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
## prepare Hoff Phantom

Hoff2 <- -!Hoff
Hoff2[is.na(Hoff2)] <- -2
mask <- !is.na(Hoff)
in.mask.rows <- (1:nrow(Hoff2))[X = apply(Hoff2 ==0, MARGIN = 1, FUN = any)]
in.mask.cols <- (1:ncol(Hoff2))[X = apply(Hoff2 ==0, MARGIN = 2, FUN = any)]


## run AR-FAST (Robust)
mask <- array(mask,dim=dim(tmap))
ff <- FASTfMRI(spm = tmap,method = "AR",mask = mask,verbose=TRUE,two.sided=FALSE)
tmap1 <- tmap
tmax <- max(abs(tmap))
tmap[ff$ActMap==0] <- NA


par(mfrow=c(3,5),mar=rep(0,4))
for(i in 1:15){
image(Hoff[in.mask.rows,in.mask.cols],col=grey(1-c(0,0.05,0.025)),axes=FALSE)
image(tmap[in.mask.rows,in.mask.cols,i],col=burd(16),zlim=c(-tmax,tmax),axes=FALSE,add=TRUE)
contour(Hoff2[in.mask.rows, in.mask.cols], drawlabels=FALSE, method = "simple", nlevels = 1, add = TRUE, lwd = 0.5)
}

## run AM-FAST (Likelihood)

ff <- FASTfMRI(spm = tmap1,method = "AM",mask = mask, verbose=TRUE,two.sided=FALSE)

tmap1[ff$ActMap==0] <- NA

par(mfrow=c(3,5),mar=rep(0,4))
for(i in 1:15){
image(Hoff[in.mask.rows,in.mask.cols],col=grey(1-c(0,0.05,0.025)),axes=FALSE)
image(tmap1[in.mask.rows,in.mask.cols,i],col=burd(16),zlim=c(-tmax,tmax),axes=FALSE,add=TRUE)
contour(Hoff2[in.mask.rows, in.mask.cols], drawlabels=FALSE, method = "simple", nlevels = 1, add = TRUE, lwd = 0.5)
}


}

