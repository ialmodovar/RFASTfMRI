\name{Jaccard Index similarity coefficient}
\alias{jaccard.index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Jaccard Index similarity coefficient
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
jaccard.index(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of 0 and 1.
%%     ~~Describe \code{x} here~~
}
  \item{y}{ vector of 0 and 1.
%%     ~~Describe \code{y} here~~
}
}
\details{
x and y need to be of the same size. NA will be automatically remove. Value close to 1 means excellent similarity while value close to 0 vector are different
}

\references{
Almodóvar-Rivera, I., & Maitra, R. (2017). FAST adaptive smoothing and thresholding for improved activation detection in low-signal fMRI. arXiv preprint arXiv:1702.00111.

Maitra, R. (2010). A re-defined and generalized percent-overlap-of-activation measure for studies of fMRI reproducibility and its use in identifying outlier activation maps. Neuroimage, 50(1), 124-135.

Jaccard, P. (1901). Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bull Soc Vaudoise Sci Nat, 37, 547-579.
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

