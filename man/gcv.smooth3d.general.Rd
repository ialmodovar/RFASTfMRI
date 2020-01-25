\name{Robust Smoothing 3D}
\alias{gcv.smooth3d.general}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust Smoothing 3D
}
\description{
Perform robust smoothing of Garcia (2010).
}
\usage{
gcv.smooth3d.general(y, initval)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
3D array
}
  \item{initval}{ Range of initial value
}
}

\value{
im.smooth: smoothed array

par: values of the estimate smooth
}
\references{
Garcia, D. (2010). Robust smoothing of gridded data in one and higher dimensions with missing values. Computational statistics & data analysis, 54(4), 1167-1178.
}
\author{
Israel Almodovar-Rivera and Ranjan Maitra.
}

