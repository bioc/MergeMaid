\name{intcorDens}
\alias{intcorDens}
\title{plot of density functions of integrative correlations}

\description{Given a mergeExprSet, this function calculates and plots the density function of the integrative correlations, as well as densities for two "null distributions" obtained by randomly permuting sample IDs.}

\usage{
   intcorDens(x,method= c("pearson", "kendall", "spearman"),...)
}

\arguments{
   \item{x}{Object of class mergeExprSet.}
   \item{method}{Method used to calculate correlation coefficient.}
   \item{...}{graphical parameters to be passed to plot.}
}

\value{The value is null.  Returns a plot.}
\details{Integrative correlation coefficients are calcualted as follows.  The first step is to identify the n genes common to all studies.
  Within each study, we calculate the correlation coefficient between
  gene g, and every other common gene.  This gives a vector of length
  n-1.  For a pair of studies, S1 and S2, we calculate the correlation
  of correlations for gene g.  When there are more than 2 studies under
  consideration, all pairwise correlation of correlations are calculated
  and averaged.  }

\seealso{\code{\link{mergeCor-class}},\code{\link{intCor}},\code{\link{modelOutcome}}}

\examples{
 
  data(mergeData)
  merged  <-mergeExprs(sample1,sample2,sample3)
  intcorDens(merged)

  intcorDens(merged,cex.legend=1.5)

  intcorDens(merged,lty=2)
  
}

\keyword{hplot}
