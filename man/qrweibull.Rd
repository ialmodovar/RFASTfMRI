\name{Reverse Weibull distribution}
\alias{qrweibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reverse Weibull distribution
}
\description{
Quantile function for the reverse (or negative) Weibull distribution with location, scale and shape parameters}
\usage{
qrweibull(p, mu = 0, beta = 1, alpha = 1, lower.tail = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
vector of probabilities
}
  \item{mu}{
Location parameter (can be given as a vector).
}
  \item{beta}{
Scale parameter (can be given as a vector).
}
  \item{alpha}{
Shape parameter (can be given as a vector).
}
  \item{lower.tail}{
Logical; if TRUE (default), probabilities are P(X<=x), otherwise P(X>x).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
qrweibull gives the quantile function.

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{qgumbel}}
}
