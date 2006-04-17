\name{mergeExpressionSet}
\docType{class}
\alias{mergeExpressionSet-class}
\alias{mergeExpressionSet}

\alias{geneStudy}
\alias{geneStudy<-}
\alias{subsetmES}
\alias{geneStudy,mergeExpressionSet-method}
\alias{exprs,mergeExpressionSet-method}
\alias{names,mergeExpressionSet-method}
\alias{notes,mergeExpressionSet-method}
\alias{phenoData,mergeExpressionSet-method}
\alias{length,mergeExpressionSet-method}
\alias{summary,mergeExpressionSet-method}
\alias{geneNames,mergeExpressionSet-method}

\alias{exprs<-,mergeExpressionSet-method}
\alias{geneNames<-,mergeExpressionSet-method}
\alias{names<-,mergeExpressionSet-method}
\alias{notes<-,mergeExpressionSet-method}
\alias{phenoData<-,mergeExpressionSet-method}
\alias{geneStudy<-,mergeExpressionSet-method}
\alias{[,mergeExpressionSet-method}

\alias{plot,mergeExpressionSet-method}
\alias{intcorDens,mergeExpressionSet-method}
\alias{intCor,mergeExpressionSet-method}
\alias{intersection,mergeExpressionSet-method}
\alias{modelOutcome,mergeExpressionSet-method}

\alias{isna}
\alias{maxintcor}


\title{Class mergeExpressionSet, a class for merged microarray data, and methods for
processing them}

\description{This is class representation for merged Microarray Data.}

\section{Slots}{
   We assume there are K studies, representing a total of M unique genes.
   \describe{
      \item{data}{A list of ExpressionSet objects.}		      
      \item{geneStudy}{Binary incidence matrix with M rows and K columns. Each
       column represents a study, and each row represents a gene. If study "s" contains gene "g", then geneStudy[g,s]=1, otherwise geneStudy[g,s]=0.
      }
      \item{notes}{Object of class "character" This slot is available for storage of descriptive information.}  
   }
}


\section{Methods}{
   Derived from ExpressionSet:
   \describe{
      \item{exprs(mergeExpressionSet)}{ An accessor function for the data slot.} 
      \item{exprs<- (mergeExpressionSet)}{A replace function for the data slot}
      \item{notes (mergeExpressionSet)}{An accessor function for the notes slot.} 
      \item{notes<- (mergeExpressionSet)}{A replace function for the notes slot.}
      \item{geneNames (mergeExpressionSet)}{Accessor function for union of gene ids in all studies.} 
      \item{geneNames<- (mergeExpressionSet)}{A replace function for gene ids.} 
 }
   Class-specific methods:
   \describe{
      \item{geneStudy (mergeExpressionSet)}{Accessor function for the geneStudy slot.}
      \item{phenoData (mergeExpressionSet)}{Accessor function for phenodata in ExpressionSet's. Returns a list, one phenodata matrix per study.} 
      \item{intersection (mergeExpressionSet)}{Represent data for genes common to all studies as a single ExpressionSet object.}
      \item{modelOutcome (mergeExpressionSet)}{Calculate regression coefficients for each study/gene.}
      \item{intCor (mergeExpressionSet)}{Calculate the integrative correlation coefficients for mergeExpressionSet data.}
      \item{plot (mergeExpressionSet)}{Draw scatterplots to compare integrative correlations for genes, here we use the approximate method to calculate the integrative correlation.}
      \item{intcorDens (mergeExpressionSet)}{Plot the distribution of  the integrative correlation coefficients and the null distribution obtained by permutation here we use the approximate method to calculate the integrative correlation.}
  }
  Standard generic methods:
  \describe{
      \item{length (mergeExpressionSet)}{Function returning the number of studies in the mergeExpressionSet.} 
      \item{names (mergeExpressionSet)}{Function returning study names.} 
      \item{names<- (mergeExpressionSet)}{A replace function for study names.}
      \item{[ (mergeExpressionSet)}{ A subset operator. Returns a mergeExpressionSet containing a subset of the studies.  A mergeExpressionSet with only one study is returned as a single ExpressionSet.}
      \item{summary (mergeExpressionSet)}{Obtain the basic information for 'mergeExpressionSet'.}

  }
}



\details{The mergeExpressionSet class is conceived as an extension of the ExpressionSet class provided in Biobase for the storage of expression array data.  A mergeExpressionSet object is primarily a list of ExpressionSet objects, along with an incidence matrix indicating which genes appear in which studies.  A mergeExpressionSet object with a single study reverts to the ExpressionSet class. A number of accessor functions are defiined for this class, as well as a few convenient analysis and plotting functions. }

\seealso{
  \code{\link{mergeExprs}},\code{\link{intCor}},\code{\link{modelOutcome}},\code{\link{intcorDens}}}

\examples{
  if(require(Biobase) & require(MASS)){
  data(mergeData)
  merged  <-mergeExprs(sample1,sample2,sample3)

  merged[1:2]
  i<-c(1,3)
  merged[i]

  exprs(merged)

  names(merged)<-c("study1","study2","study3")

  length(merged)

  summary(merged)

  plot(merged)

  plot(merged[1:2])

  intcorDens(merged)

  inter  <- intersection(merged)
  }
}

\keyword{classes}
