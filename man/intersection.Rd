\name{intersection}
\alias{intersection}
\title{eSet with all common genes in a mergeESet}

\description{Given a mergeESet, this function returns a single eSet.  Only genes common to all studies are included.  Expression data for all studies sits side by side in the 'exprs' slot. The 'notes' slot is used for inormation about the study identity of each sample.}

\usage{
   intersection(x)
}

\arguments{
   \item{x}{Object of class mergeESet.}
   
}

\value{Returns an object of class eSet}


\seealso{  \code{\link{mergeESet-class}}}


\examples{
  
  data(mergeData)
  merged  <-mergeExprs(sample1,sample2,sample3)

  inter  <- intersection(merged)
  note  <- notes(inter)
  
}

\keyword{manip}
