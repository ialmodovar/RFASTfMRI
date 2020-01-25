\name{Gauss.smooth}
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
  \item{scaled}{If TRUE (default) scaled the smooth }
}

\value{
 Return smoothed map using the provided FWHM.
}
\references{
Almodovar-Rivera, I., & Maitra, R. (2019). FAST adaptive smoothing and thresholding for improved activation detection in low-signal fMRI. IEEE transactions on medical imaging, 38(12), 2821-2828.
}
\author{
Israel Almodovar-Rivera and Ranjan Maitra.
}

\examples{
\dontrun{
data(tmap)
##before smoothing
t1 <- tmap[,,1]
image(t1,main="Before smooting")
sm <- Gauss.smooth(tstat = t1, fwhm=c(3,3))
image(sm,main="Smoothing with FWHM (3,3)")
}
}