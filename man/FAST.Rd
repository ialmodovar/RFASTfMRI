\name{Fast Adaptive Smoothing and Thresholding}
\alias{FAST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Adaptive Smoothing and Thresholding
}
\description{
Compute FAST for activation-detection in fMRI studies.By default this one perform one sided tes
}
\usage{
FAST(spm, method = "robust", mask = NULL, alpha = 0.05, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spm}{ Statistical Parametric Map. Array of 3D.
}
  \item{method}{If method = "robust" will perform robust smoothing of Garcia (2010). If method = "likelihood" will perform maximum estimate of the optimal FWHM to perform smoothing.
}
  \item{mask}{
Region of interest. If not mask is provided a estimate one will be obtained
}
  \item{alpha}{Significance level, default 0.05
}
  \item{verbose}{ If TRUE will print every calculation step. Default is FALSE
}
}
\details{
Compute FAST for activation detection in fMRI studies. By default this one perform one sided test.
}
\value{
FWHM: a 3D array containing the value of the estimated FWHM at each iteration.

ActMap: a 3D array of 0 and 1. If a voxel is determined to activated it is classified as 1.

SPM: initial statistical parametric map.

SmoothSPM: Final smoothed map.

JaccardIndex: Jaccard index at each iteration

Rho: Estimate of the spatial correlation at each step

An: Normalizing sequence of the extreme value theory.

Bn: Normalizing sequence of the extreme value theory.

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
