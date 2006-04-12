\name{intCor}
\alias{intCor}
\title{Correlation of Correlations}

\description{Given a mergeESet, this function calculates the study specific
  correlation matrices, and, for each gene, the correlation of
  correlations. }

\usage{
   intCor(x,method= c("pearson", "spearman"),exact,...)
}

\arguments{
   \item{x}{Object of class mergeESet.}
   \item{method}{Method used to calculate correlation coefficient. If exact is TRUE, the available methods to use is "spearman" and "pearson"; If exact is FALSE, the available methods to use is "pearson".}
  \item{exact}{If exact is TRUE, we use the standard method the calculate  the integrative correlation; If exact is FALSE, we use the approximate method the calculate.}
\item{...}{Not implemented at this time}
}

\value{The output is an object of class  mergeCor.}
\details{Integrative correlation coefficients are calcualted as follows.  The first step is to identify the n genes common to all studies.
Within each study, we calculate the correlation coefficient between gene g, and every other common gene.  This gives a vector of length n-1.  For a pair of studies, S1 and S2, we calculate the correlation of correlations for gene g.  When there are more than 2 studies under consideration, all pairwise correlation of correlations are calculated and averaged.  }


\seealso{
  \code{\link{mergeCor-class}},\code{\link{intcorDens}}}

\examples{

if(require(Biobase) & require(MASS)){
  data(mergeData)
  merged  <-mergeExprs(sample1,sample2,sample3)
  corcor  <-intCor(merged,method="spearman")
  
  plot(merged)
  hist(corcor)

  corcor  <-intCor(merged,method="pearson",exact=FALSE)
  corcor  <-intCor(merged[1:2])
  corcor  <-intCor(merged,exact=TRUE)

  vv<-c(1,3)
  corcor1  <-intCor(merged[vv])
  plot(merged,xlab="study A",ylab="study B",main="CORRELATION OF CORRELATION",col=3,pch=4)
  hist(corcor1,xlab="CORRELATION OF CORRELATION")
 }  
}

\keyword{univar}
