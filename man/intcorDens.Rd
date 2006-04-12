\name{intcorDens}
\alias{intcorDens}
\title{plot of density functions of integrative correlations}

\description{Given a mergeESet, this function calculates and plots the density function of the approximate integrative correlations, as well as densities for the "null distributions" obtained by randomly permuting sample IDs.}

\usage{
   intcorDens(x,method,...)
}

\arguments{
   \item{x}{Object of class mergeESet.}
   \item{method}{The available method to use is "pearson".}
   \item{...}{Graphical parameters to be passed to plot.}
}

\value{The value is null.  Returns a plot.}
\details{Here we use the approximate method to calculate the integrative correlation.  }

\seealso{\code{\link{mergeESet-class}},\code{\link{intCor}},\code{\link{modelOutcome}}}

\examples{
  if(require(Biobase) & require(MASS)){
  data(mergeData)
  merged  <-mergeExprs(sample1,sample2,sample3)

  intcorDens(merged)

  intcorDens(merged,cex.legend=1.5)

  intcorDens(merged,lty=2)
  }
}

\keyword{hplot}
