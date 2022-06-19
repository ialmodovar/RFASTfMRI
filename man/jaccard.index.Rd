\name{jaccard.index}
\alias{jaccard.index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Jaccard Index similarity coefficient
}
\description{
Compute Jaccard index to determine similarity between activation
maps. 
}
\usage{
jaccard.index(x, y, value.if.both.nulls = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of 0's (inactives) and 1's (actives).
%%     ~~Describe \code{x} here~~
}
  \item{y}{ vector of 0's (inactives) and 1's (actives).
%%     ~~Describe \code{y} here~~
}
  \item{value.if.both.nulls}{value to be returned if both x and y are all 0. By default, this is set to 1 but can be set by the user to a desired value, for instance, when it is known that x and y are incorrectly all 0.
}
\details{
x and y need to be of the same size. NAs will be automatically removed. A value close to 1 means excellent similarity between the actives in x and those in y while values close to zero indicate very little similarity.
}
}

\references{
Almodovar-Rivera, I., & Maitra, R. (2019). FAST adaptive smoothing and thresholding for improved activation detection in low-signal fMRI. IEEE Transactions on Medical Imaging, 38(12), 2821-2828.

Maitra, R. (2010). A re-defined and generalized percent-overlap-of-activation measure for studies of fMRI reproducibility and its use in identifying outlier activation maps. Neuroimage, 50(1), 124-135.

Jaccard, P. (1901). Etude comparative de la distribution florale dans une portion des Alpes et des Jura. Bull Soc Vaudoise Sci Nat, 37, 547-579.
}
\author{
Israel Almodovar-Rivera and Ranjan Maitra.
}

\examples{
\dontrun{
set.seed(787)
x <- sample(x = 0:1,size = 50,replace=TRUE)
y <- sample(x = 0:1,size = 50,replace=TRUE)
jaccard.index(x = x,y = y)
}
}
