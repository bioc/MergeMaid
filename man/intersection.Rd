\name{intersection}
\alias{intersection}
\title{exprSet with all common genes in a mergeExprSet}

\description{Given a mergeExprSet, this function returns a single exprSet.  Only genes common to all studies are included.  Expression data for all studies sits side by side in the 'exprs' slot. The 'notes' slot is used for inormation about the study identity of each sample.}

\usage{
   intersection(x)
}

\arguments{
   \item{x}{Object of class mergeExprSet.}
   
}

\value{Returns an object of class exprSet}


\seealso{  \code{\link{mergeExprSet-class}}}


\examples{

  data(mergeData)
  merged  <-mergeExprs(sample1,sample2,sample3)

  inter  <- intersection(merged)
  note  <- notes(inter)
  
}

\keyword{manip}
