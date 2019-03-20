\name{FAST}
\alias{FAST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Adaptive Smoothing and Thresholding for low signal fMRI
}
\description{
Compute FAST for activation-detection in fMRI studies.
}
\usage{
FAST(spm, method = "robust", mask = NULL, alpha = 0.05, gcv.init=NULL,
verbose = FALSE,all = FALSE,stopping=TRUE, K = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spm}{ Statistical Parametric Map. Array of size 2D or 3D.
}
  \item{method}{If method = "robust" will perform robust smoothing of Garcia (2010). If method = "likelihood" will perform maximum estimate of the optimal FWHM to perform smoothing.
}
  \item{mask}{Region of interest. If not mask is provided a estimate one
    will be obtained.}
   \item{alpha}{Significance level, default 0.05.}
  \item{gcv.init}{ If NULL, will choose initial values for optim() from
    a range of values than maximized the gcv score function. If provided the length of the FWHM must be from the same as dimension as the image. }
  \item{verbose}{ If TRUE, will print every calculation step. Default is
    FALSE.}
   \item{all}{Return everything FAST compute. Default is activated
     map, smoothed spm and original spm. }
  \item{stopping}{If TRUE algorithm will use stopping criterion
    suggested in  Almodóvar-Rivera and Maitra (2018).}
  \item{K}{Number of desired iterations only works if stopping is FALSE.}
}
\details{
Compute FAST for activation detection in fMRI studies. By default perform one sided test.
}
\value{
FWHM: a 2D or 3D array containing the value of the estimated FWHM at each iteration.

ActMap: a 2D or 3D array of 0 and 1. If a voxel is determined to activated it is classified as 1.

SPM: initial statistical parametric map.

SmoothSPM: Final smoothed t-map.

JaccardIndex: Jaccard index between activation maps at each iteration.

Rho: Estimate of the spatial correlation at each step.

An: Normalizing sequence of the extreme value theory.

Bn: Normalizing sequence of the extreme value theory. 

}
\references{
Almodóvar-Rivera, I., & Maitra, R. (2018). FAST adaptive smoothing and thresholding for improved activation detection in low-signal fMRI. arXiv preprint arXiv:1702.00111.

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
