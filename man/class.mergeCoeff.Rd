\name{mergeCoeff}
\alias{mergeCoeff-class}
\alias{mergeCoeff}
\alias{check.length}
\alias{coeff}
\alias{coeff<-}
\alias{stdcoeff}
\alias{stdcoeff<-}
\alias{zscore}
\alias{zscore<-}
\alias{coeff,mergeCoeff-method}
\alias{stdcoeff,mergeCoeff-method}
\alias{zscore,mergeCoeff-method}
\alias{coeff<-,mergeCoeff-method}
\alias{stdcoeff<-,mergeCoeff-method}
\alias{zscore<-,mergeCoeff-method}

\alias{mergeCoeff}
\title{Class mergeCoeff, a class for storing regression coefficients.}

\description{This is the class representation for output from regression coefficient calculations} 

\section{Slots}{
   \describe{
      \item{coeff}{Object of class matrix, A matrix storing default
     coefficients. }
      \item{coeff.std}{Object of class matrix, A matrix storing standardized
     coefficients.}
     \item{zscore}{Object of class matrix, A matrix storing zscores.}	
 }
}


\section{Methods}{
   Class-specific methods:
   \describe{
      \item{coeff (mergeCoeff)}{Accessor function for the coeff slot.} 
      \item{coeff<- (mergeCoeff)}{Replacement function for the coeff slot.}
      \item{stdcoeff (mergeCoeff)}{Accessor function for the coeff.std slot.} 
      \item{stdcoeff<- (mergeCoeff)}{Replacement function for the coeff.std slot.}
      \item{zscore (mergeCoeff)}{Accessor function for the zscore slot.} 
      \item{zscore<- (mergeCoeff)}{Replacement function for the zscore slot.}
  }
  Standard generic methods:
  \describe{
      \item{plot (list)}{This method is not formally defined for mergeCoeff objects but for a matrix.  This function would typically be called with the following syntax, plot(coeff(mergeCoeff)).The result is pairwise scatterplots of the columns of the selected matrix.  If there are two studies, this is a single scatterplot.}
  }
}


\seealso{ \code{\link{mergeExprs}},\code{\link{modelOutcome}}, \code{\link{mergeExpressionSet-class}}}


\examples{
  if(require(Biobase) & require(MASS) & require(survival)){
  
  data(mergeData)
  merged  <- mergeExprs(sample1,sample2,sample3)

  log.coeff  <- modelOutcome(merged,outcome=c(1,1,1),method="logistic")
  plot(coeff(log.coeff))
  plot(stdcoeff(log.coeff),pch=4,labels=c("study A","study B","study C"),col=3)

  linear.coeff  <- modelOutcome(merged[1:2],outcome=c(3,3),method="linear")
  plot(zscore(linear.coeff))
  plot(zscore(linear.coeff),xlab="study A",ylab="study B",col=2)
  }
}

\keyword{classes}
