\name{Dataset}
\alias{Hoff}
\alias{tmap}
\alias{HoffActiv}
\docType{data}
\title{Dataset for Hoff Phantom}
\description{
  Small datasets to test and demonstrate \pkg{RFASTfMRI}.
}
\usage{
data(Hoff)
data(HoffActiv)
data(tmap)
}
\format{
  \code{Hoff}, \code{HoffActiv} are a 2D dimensional array of size 128  x 128.
  
  \code{tmap} is a 3D dimensional array of size 128 x 128 x 15.

 \code{HoffActiv} is 128 x 128 array. .
}
\details{
  \code{Hoff} is 128 x 128 array. Voxels that belong to the region-of-interest are classify as 0 and 1, for a total 3485. Outside this region is NA.
  
  \code{HoffActiv} is 128 x 128 array. Voxels that are activated are
  classify as 2 for a total 138. Outside this region is NA.

  \code{tmap} a 128 x 128 x 15 array. The third coordinate are independent replicates, a total of 15 maps were simulate.
}

\author{
Israel A. Almodovar-Rivera and Ranjan Maitra
}
\keyword{datasets}