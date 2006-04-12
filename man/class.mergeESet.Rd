\name{mergeESet}
\docType{class}
\alias{mergeESet-class}
\alias{mergeESet}

\alias{geneStudy}
\alias{geneStudy<-}
\alias{subsetmES}
\alias{geneStudy,mergeESet-method}
\alias{exprs,mergeESet-method}
\alias{names,mergeESet-method}
\alias{notes,mergeESet-method}
\alias{phenoData,mergeESet-method}
\alias{length,mergeESet-method}
\alias{summary,mergeESet-method}
\alias{geneNames,mergeESet-method}

\alias{exprs<-,mergeESet-method}
\alias{geneNames<-,mergeESet-method}
\alias{names<-,mergeESet-method}
\alias{notes<-,mergeESet-method}
\alias{phenoData<-,mergeESet-method}
\alias{geneStudy<-,mergeESet-method}
\alias{[,mergeESet-method}

\alias{plot,mergeESet-method}
\alias{intcorDens,mergeESet-method}
\alias{intCor,mergeESet-method}
\alias{intersection,mergeESet-method}
\alias{modelOutcome,mergeESet-method}

\alias{isna}
\alias{maxintcor}


\title{Class mergeESet, a class for merged microarray data, and methods for
processing them}

\description{This is class representation for merged Microarray Data.}

\section{Slots}{
   We assume there are K studies, representing a total of M unique genes.
   \describe{
      \item{data}{A list of eSet objects.}		      
      \item{geneStudy}{Binary incidence matrix with M rows and K columns. Each
       column represents a study, and each row represents a gene. If study "s" contains gene "g", then geneStudy[g,s]=1, otherwise geneStudy[g,s]=0.
      }
      \item{notes}{Object of class "character" This slot is available for storage of descriptive information.}  
   }
}


\section{Methods}{
   Derived from \code{\link{eSet}}:
   \describe{
      \item{exprs(mergeESet)}{ An accessor function for the data slot.} 
      \item{exprs<- (mergeESet)}{A replace function for the data slot}
      \item{notes (mergeESet)}{An accessor function for the notes slot.} 
      \item{notes<- (mergeESet)}{A replace function for the notes slot.}
      \item{geneNames (mergeESet)}{Accessor function for union of gene ids in all studies.} 
      \item{geneNames<- (mergeESet)}{A replace function for gene ids.} 
 }
   Class-specific methods:
   \describe{
      \item{geneStudy (mergeESet)}{Accessor function for the geneStudy slot.}
      \item{phenoData (mergeESet)}{Accessor function for phenodata in eSet's. Returns a list, one phenodata matrix per study.} 
      \item{intersection (mergeESet)}{Represent data for genes common to all studies as a single eSet object.}
      \item{modelOutcome (mergeESet)}{Calculate regression coefficients for each study/gene.}
      \item{intCor (mergeESet)}{Calculate the integrative correlation coefficients for mergeESet data.}
      \item{plot (mergeESet)}{Draw scatterplots to compare integrative correlations for genes, here we use the approximate method to calculate the integrative correlation.}
      \item{intcorDens (mergeESet)}{Plot the distribution of  the integrative correlation coefficients and the null distribution obtained by permutation here we use the approximate method to calculate the integrative correlation.}
  }
  Standard generic methods:
  \describe{
      \item{length (mergeESet)}{Function returning the number of studies in the mergeESet.} 
      \item{names (mergeESet)}{Function returning study names.} 
      \item{names<- (mergeESet)}{A replace function for study names.}
      \item{[ (mergeESet)}{ A subset operator. Returns a mergeESet containing a subset of the studies.  A mergeESet with only one study is returned as a single eSet.}
      \item{summary (mergeESet)}{Obtain the basic information for 'mergeESet'.}

  }
}



\details{The mergeESet class is conceived as an extension of the eSet class provided in Biobase for the storage of expression array data.  A mergeESet object is primarily a list of eSet objects, along with an incidence matrix indicating which genes appear in which studies.  A mergeESet object with a single study reverts to the eSet class. A number of accessor functions are defiined for this class, as well as a few convenient analysis and plotting functions. }

\seealso{
  \code{\link{mergeExprs}},\code{\link{intCor}},\code{\link{modelOutcome}},\code{\link{intcorDens}},eSet}

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
  note  <- notes(inter)
  }
}

\keyword{classes}
