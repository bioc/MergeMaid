\name{intersection}
\alias{intersection}
\title{ExpressionSet with all common genes in a mergeExpressionSet}

\description{Given a mergeExpressionSet, this function returns a single ExpressionSet.  Only genes common to all studies are included.  Expression data for all studies sits side by side in the 'exprs' slot. The 'notes' slot is used for inormation about the study identity of each sample.}

\usage{
   intersection(x)
}

\arguments{
   \item{x}{Object of class mergeExpressionSet.}
   
}

\value{Returns an object of class ExpressionSet}


\seealso{  \code{\link{mergeExpressionSet-class}}}


\examples{
  
  data(mergeData)
  merged  <-mergeExprs(sample1,sample2,sample3)

  inter  <- intersection(merged)
 
}

\keyword{manip}
