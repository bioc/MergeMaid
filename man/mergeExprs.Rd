\name{mergeExprs}
\alias{mergeExprs}
\alias{AverageDuplicates}
\alias{check}
\alias{mergeget}

\title{Merge gene expression data sets}

\description{Merges gene expression data from different studies.}

\usage{
   mergeExprs(...)
}

\arguments{
   \item{...}{Input objects can be any combination of mergeExprSet, exprSet, matrix or a list. 
A list should have the following slots:
expression matrix,
pheno data matrix, 
gene names vector,
notes. 
The order of the four slots is fixed.
A matrix should have genes ids as its row names, as should the exprs slot of an exprSet. 
Since merging depends on geneids, these conventions are essential.}
   
}

\value{The output is a mergeExprSet.}

\details{The mergeExprSet object is the standard input for all functions in the MergeMaid package. Use the mergeExprs function when creating mergeExprSet objects to ensure that all necessary information is available for further analysis. }

\seealso{\code{\link{mergeExprSet-class}}}

\examples{
 
  data(mergeData)
  merged  <-mergeExprs(sample1,sample2,sample3)

  rr<-rnorm(200*22,0,1)
  mm<-matrix(rr,200,22)
  rownames(mm)<-geneNames(sample2)
  merge.m<-mergeExprs(sample1,mm,sample2)
  intcor.m<-intCor(merge.m)
  plot(intcor.m)

  rr<-rnorm(200*50,0,1)
  mm2<-matrix(rr,200,50)
  ph.ll<-as.data.frame(rbinom(200,1,.5))
  ll<-list(mm2,ph.ll,geneNames(sample2),"list 2")
  merge.t<-mergeExprs(sample1,mm,sample2,ll)
  intcor.t<-intCor(merge.t)
  plot(intcor.t)

  merge.a<-mergeExprs(sample3,merge.m,ll)
  inter<-intersection(merge.a)
  summary(merge.a)
  
}

\keyword{manip}
