\name{qrweibull}
\alias{qrweibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantile of a Reverse Weibull distribution
}
\description{
Quantile function for the reverse (or negative) Weibull distribution with location, scale and shape parameters}
\usage{
qrweibull(p, mu = 0, beta = 1, tau = 1, lower.tail = TRUE)
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
  \item{tau}{
Shape parameter (can be given as a vector).
}
  \item{lower.tail}{
Logical; if TRUE (default), probabilities are P(X<=x), otherwise P(X>x).
}
}
\details{
Compute the quantile of a Reverse Weibull distribution. Accept general
location parameter mu, and shapes parameters beta and tau.
}
\value{
qrweibull gives the quantile function.

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Israel Almodóvar-Rivera and Ranjan Maitra.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{qgumbel}}
}
