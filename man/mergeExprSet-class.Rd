\name{mergeExprSet-class}
\docType{class}
\alias{mergeExprSet-class}
\alias{mergeExprSet}

\alias{geneStudy}
\alias{geneStudy<-}
\alias{subsetmES}
\alias{geneStudy,mergeExprSet-method}
\alias{exprs,mergeExprSet-method}
\alias{names,mergeExprSet-method}
\alias{notes,mergeExprSet-method}
\alias{phenoData,mergeExprSet-method}
\alias{length,mergeExprSet-method}
\alias{summary,mergeExprSet-method}
\alias{geneNames,mergeExprSet-method}

\alias{exprs<-,mergeExprSet-method}
\alias{geneNames<-,mergeExprSet-method}
\alias{names<-,mergeExprSet-method}
\alias{notes<-,mergeExprSet-method}
\alias{phenoData<-,mergeExprSet-method}
\alias{geneStudy<-,mergeExprSet-method}
\alias{[,mergeExprSet-method}

\alias{intcorDens,mergeExprSet-method}
\alias{intCor,mergeExprSet-method}
\alias{intersection,mergeExprSet-method}
\alias{modelOutcome,mergeExprSet-method}

\title{Class mergeExprSet, a class for merged microarray data, and methods for
processing them}

\description{This is class representation for merged Microarray Data.}

\creatingobjects{new('mergeExprSet',
                      data           = ..., # Object of class list,
		      geneStudy      = ..., # Object of class matrix,
		      notes          = ..., # Object of class vector.)}

\slots{We assume there are K studies, representing a total of M unique genes.
			
   \item{data}{A list of exprSet objects.}		      
   \item{geneStudy}{Binary incidence matrix with M rows and K columns. Each
     column represents a study, and each row represents a gene. If study "s" contains gene "g", then geneStudy[g,s]=1, otherwise geneStudy[g,s]=0.
    }
   \item{notes}{Object of class "character" This slot is available for storage of descriptive information.}  
}

\methods{
  
  \item{[ (mergeExprSet)}{ A subset operator. Returns a mergeExprSet containing a subset of the studies.  A mergeExprSet with only one study is returned as a single exprSet.}
  \item{exprs(mergeExprSet)}{ An accessor function for the data slot.} 

  \item{exprs<- (mergeExprSet)}{A replace function for the data slot}
  
  \item{notes (mergeExprSet)}{An accessor function for the notes slot.} 

  \item{notes<- (mergeExprSet)}{A replace function for the notes slot.}

 \item{geneStudy (mergeExprSet)}{Accessor function for the geneStudy slot.}

  \item{length (mergeExprSet)}{Function returning the number of studies in the mergeExprSet.} 

  \item{geneNames (mergeExprSet)}{Accessor function for union of gene ids in all studies.} 

  \item{geneNames<- (mergeExprSet)}{A replace function for gene ids.} 

  \item{names (mergeExprSet)}{Function returning study names.} 

  \item{names<- (mergeExprSet)}{A replace function for study names.}

  \item{phenoData (mergeExprSet)}{Accessor function for phenodata in exprSet's. Returns a list, one phenodata matrix per study.} 
 

  \item{summary (mergeExprSet)}{Obtain the basic information for 'mergeExprSet'.}

\item{intersection (mergeExprSet)}{Represent data for genes common to all studies as a single exprSet object.}

\item{modelOutcome (mergeExprSet)}{Calculate regression coefficients for each study/gene.}

   \item{intCor (mergeExprSet)}{Calculate the integrative correlation coefficients for mergeExprSet data.}

  \item{intcorDens (mergeExprSet)}{Plot the distribution of  the integrative correlation coefficients and two null distributions obtained by permutation.}
}

\Details{The mergeExprSet class is conceived as an extension of the exprSet class provided in Biobase for the storage of expression array data.  A mergeExprSet object is primarily a list of exprSet objects, along with an incidence matrix indicating which genes appear in which studies.  A mergeExprSet object with a single study reverts to the exprSet class. A number of accessor functions are defiined for this class, as well as a few convenient analysis and plotting functions. }

\seealso{
  \code{\link{mergeExprs}},\code{\link{intCor}},\code{\link{modelOutcome}},\code{\link{intcorDens}},exprSet}

\examples{
  
  data(mergeData)
  merged  <-mergeExprs(sample1,sample2,sample3)

  merged[1:2]
  i<-c(1,3)
  merged[i]

  exprs(merged)

  names(merged)<-c("study1","study2","study3")

  length(merged)

  summary(merged)

  intcorDens(merged)

  inter  <- intersection(merged)
  note  <- notes(inter)
  
}

\keyword{classes}
