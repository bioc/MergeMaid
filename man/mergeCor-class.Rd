\name{mergeCor}
\alias{mergeCor-class}
\alias{mergeCor}

\alias{cors}
\alias{pairwise.cors}
\alias{integrative.cors}
\alias{cors,mergeCor-method}
\alias{pairwise.cors,mergeCor-method}
\alias{integrative.cors,mergeCor-method}
\alias{plot,mergeCor-method}
\alias{hist,mergeCor-method}
\alias{plot,list-method}


\title{Class mergeCor, a class for storing data relevant to integrative correlation coefficients.}

\description{This is the class representation for integrative correlation coefficients.} 

\creatingobjects{new('mergeCor',
                      cors = ..., # Object of class list 
                      pairwise.cors   = ..., # Object of class matrix
		    )}

\slots{
			
  \item{cors}{Object of class list.  Each element is the correlation matrix for a single study.}
   \item{pairwise.cors}{Object of class matrix.  Each column contains  correlation of correlation score for genes in a pair of study.}	
}

\methods{
\item{cors(mergeCor)}{Accessor function for the cors slot.}
\item{pairs.cors( mergeCor)}{Accessor function for the pairwise.cors slot}
\item{integrate.cors(mergeCor)}{Accessor function, returns integrative correlation (average pairwise correlation of correlations) for each gene.}
\item{plot(mergeCor,...)}{Draw scatterplots to compare correlations of genes.  }
\item{hist(mergeCor,...)}{Draw histograms of integrative correlation.  }
}

\details{ If 'n' is the number of studies
 then for i < j <= n, the pairwise correlation of correlations for studies i and j is stored in column ((i-1)*(n-1)-(i-2)*(i-1)/2 + j-i) of the pairwise.cors slot. }

\seealso{ \code{\link{mergeCor-class}},\code{\link{intCor}},\code{\link{modelOutcome}}, \code{\link{mergeExprSet-class}}}

\examples{
  
  data(mergeData)
  merged  <- mergeExprs(sample1,sample2,sample3)
  intcor3  <-intCor(merged,method="pearson")
  plot(intcor3)
  intcor2  <-intCor(merged[1:2],method="kendall")
  plot(intcor2,pch=4,col=5)

  cors(intcor2)
  pairwise.cors(intcor3)
  integrative.cors(intcor3)
  integrative.cors(intcor2)["Hs.12101"]
  
}

\keyword{classes}
