\name{Gaussian Smooth}
\alias{Gauss.smooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gaussian Smooth 
}
\description{
Perform Gaussian smoothing in 2D and 3D array
}
\usage{
Gauss.smooth(tstat, fwhm, scaled = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tstat}{ A 2D or 3D array of t-map or (z-map).
%%     ~~Describe \code{tstat} here~~
}
  \item{fwhm}{ Vector in terms of the Full-width Half-Maximum (FWHM).
%%     ~~Describe \code{fwhm} here~~
}
  \item{scaled}{
If TRUE (default) scaled the smooth 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 Return smoothed map using the provided FWHM
}
\references{
Almodóvar-Rivera, I., & Maitra, R. (2017). FAST adaptive smoothing and thresholding for improved activation detection in low-signal fMRI. arXiv preprint arXiv:1702.00111.
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

