\name{mergeCor}
\alias{mergeCor-class}
\alias{mergeCor}

\alias{pairwise.cors}
\alias{integrative.cors}
\alias{pairwise.cors,mergeCor-method}
\alias{integrative.cors,mergeCor-method}
\alias{maxcors,mergeCor-method}
\alias{hist,mergeCor-method}
\alias{plot,list-method}
\alias{notes,mergeCor-method}

\alias{hist}
\alias{plot}
\alias{maxcors}

\title{Class mergeCor, a class for storing data relevant to integrative correlation coefficients.}

\description{This is the class representation for integrative correlation coefficients.} 

\section{Slots}{
   \describe{
      \item{pairwise.cors}{Object of class matrix.  Each column contains  correlation of correlation score for genes in a pair of study.}
      \item{maxcors}{Object of class vector. Each slot represents maximal canonical correlation (pairwise canonical correlations) for each pair of studies.}
      \item{notes}{Object of class vector. Each slot contains notes for each study.}
   }
}


\section{Methods}{
   Class-specific methods:
   \describe{
      \item{cors(mergeCor)}{Accessor function for the cors slot.}
      \item{pairs.cors( mergeCor)}{Accessor function for the pairwise.cors slot}
      \item{integrate.cors(mergeCor)}{Accessor function, returns integrative correlation (average pairwise correlation of correlations) for each gene. If adjust is TRUE, returns the integrate correlation devided by the correponding canonical correlation, and the default value is FALSE.}
      \item{maxcors(mergeCor)}{Accessor function, returns maximal canonical correlation (pairwise canonical correlations) for each pair of studies.}
  }
  Standard generic methods:
  \describe{
      \item{notes (mergeCor)}{An accessor function for the notes slot.} 
      \item{hist(mergeCor,...)}{Draw histograms of integrative correlation, here we use the approximate method to calculate the integrative correlation.  }
  }
}


\details{ If 'n' is the number of studies
 then for i < j <= n, the pairwise correlation of correlations for studies i and j is stored in column ((i-1)*(n-1)-(i-2)*(i-1)/2 + j-i) of the pairwise.cors slot. }

\seealso{ \code{\link{mergeCor-class}},\code{\link{intCor}},\code{\link{modelOutcome}}, \code{\link{mergeExpressionSet-class}}}

\examples{
  if(require(Biobase) & require(MASS)){
  data(mergeData)
  merged  <- mergeExprs(sample1,sample2,sample3)
  intcor3  <-intCor(merged,method="pearson")
  plot(merged)
  intcor2  <-intCor(merged[1:2],exact=FALSE)
  plot(merged,pch=4,col=5)

  pairwise.cors(intcor3)
  integrative.cors(intcor3)
  integrative.cors(intcor2)["Hs.12101"]
  maxcors(intcor2)
  }
}

\keyword{classes}
