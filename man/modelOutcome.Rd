\name{modelOutcome}
\alias{modelOutcome}
\title{Compare regression coefficients across studies}

\description{Given a set of merged studies, this function calculates study specific regression coefficients for each gene. }

\usage{
   modelOutcome(x,outcome,outcome2=NULL,method=c("linear","logistic","cox"),...)
}

\arguments{
   \item{x}{Object of class mergeESet.}
   \item{method}{Method specifies the model used to generate coefficients.  At this time only linear regression, logistic regression, and Cox hazard rates are implemented.}
   \item{outcome, outcome2}{The format for the outcome variable depends on the model used.  For linear regression, outcome should be a continous response variable, for logistic regression, it should be a binary response variable, and for Cox hazard rates it should be time of event.  Outcome 2 is currently used only in the calculation of hazard rates, and should be a binary variable indicating censoring status for each subject.
If outcome is a vector of length equal to number of studies, then each element represents the column in the eSet phenoData slot for that study.  If outcome is a list, then each list element should have actual outcome data for the corresponding study. }
\item{...}{Not implemented at this time}
}

\value{The output is a mergeCoeff object.}

\seealso{ \code{\link{modelOutcome}}, \code{\link{mergeCoeff-class}}}


\examples{
  if(require(Biobase) & require(MASS) & require(survival)){
  data(mergeData)
  merged  <- mergeExprs(sample1,sample2,sample3)

  log.coeff  <- modelOutcome(merged,outcome=c(1,1,1),method="logistic")
  plot(coeff(log.coeff))
  
  linear.coeff  <- modelOutcome(merged[1:2],outcome=c(3,3),method="linear")
  plot(zscore(linear.coeff),xlab="study 1",ylab="study 2")

  event1<-rbinom(100,1,.5)
  event2<-rbinom(50,1,.5) 
  event3<-rbinom(70,1,.5)

  out1<-rnorm(100,5,1)
  out2<-rnorm(50,5,1)
  out3<-rnorm(70,5,1)

  out<-list(out1,out2,out3)
  even<-list(event1,event2,event3)

  cox.coeff<-modelOutcome(merged,outcome2=even,outcome=out,method="cox")
  plot(coeff(cox.coeff))

  }
}

\keyword{models}
