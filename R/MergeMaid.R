

library(Biobase)

library(survival)

###  define class
setClass("mergeExprSet",representation(data="list",geneStudy="matrix",notes="character"))
## define two new classes for the output from corcor
## and the output from corcox, corlinear and corlogistic

setClass("mergeCor",representation(cors="list",pairwise.cors="matrix"))

setClass("mergeCoeff",representation(coeff="matrix",coeff.std="matrix",zscore="matrix",method="vector"))

### necessary generics
if(is.null(getGeneric("exprs"))) setGeneric("exprs",function(object) standardGeneric("exprs"))
if(is.null(getGeneric("geneNames"))) setGeneric("geneNames",function(object) standardGeneric("geneNames"))
if(is.null(getGeneric("names"))) setGeneric("names",function(x) standardGeneric("names"))
if(is.null(getGeneric("notes"))) setGeneric("notes",function(object) standardGeneric("notes"))
if(is.null(getGeneric("phenoData"))) setGeneric("phenoData",function(object) standardGeneric("phenoData"))
if(is.null(getGeneric("length"))) setGeneric("length",function(x) standardGeneric("length"))
if(is.null(getGeneric("summary"))) setGeneric("summary",function(object,...) standardGeneric("summary"))
if(is.null(getGeneric("geneStudy"))) setGeneric("geneStudy",function(x) standardGeneric("geneStudy"))
if(is.null(getGeneric("coeff"))) setGeneric("coeff",function(x) standardGeneric("coeff"))
if(is.null(getGeneric("stdcoeff"))) setGeneric("stdcoeff",function(x) standardGeneric("stdcoeff"))
if(is.null(getGeneric("zscore"))) setGeneric("zscore",function(x) standardGeneric("zscore"))
if(is.null(getGeneric("cors"))) setGeneric("cors",function(x) standardGeneric("cors"))
if(is.null(getGeneric("pairwise.cors"))) setGeneric("pairwise.cors",function(x) standardGeneric("pairwise.cors"))
if(is.null(getGeneric("integrative.cors"))) setGeneric("integrative.cors",function(x,...) standardGeneric("integrative.cors"))
if(is.null(getGeneric("modelOutcome"))) setGeneric("modelOutcome",function(x,outcome,outcome2=NULL,method=c("linear","logistic","cox"),...)  standardGeneric("modelOutcome"))
if(is.null(getGeneric("intCor"))) setGeneric("intCor",function(x,method=c("pearson","kendall","spearman"),...) standardGeneric("intCor"))
if(is.null(getGeneric("plot"))) setGeneric("plot",function(x,y,...) standardGeneric("plot"))
if(is.null(getGeneric("hist"))) setGeneric("hist",function(x,...)
standardGeneric("hist"))
if(is.null(getGeneric("intcorDens"))) setGeneric("intcorDens",function(x,method=c("pearson","kendall","spearman"),...) standardGeneric("intcorDens"))
if(is.null(getGeneric("intersection"))) setGeneric("intersection",function(x) standardGeneric("intersection"))


if(is.null(getGeneric("exprs<-"))) setGeneric("exprs<-",function(object,value) standardGeneric("exprs<-"))
if(is.null(getGeneric("geneNames<-"))) setGeneric("geneNames<-",function(object, value) standardGeneric("geneNames<-"))
if(is.null(getGeneric("names<-"))) setGeneric("names<-",function(x,value) standardGeneric("names<-"))
if(is.null(getGeneric("notes<-"))) setGeneric("notes<-",function(object,value) standardGeneric("notes<-"))
if(is.null(getGeneric("phenoData<-"))) setGeneric("phenoData<-",function(object,value) standardGeneric("phenoData<-"))
if(is.null(getGeneric("geneStudy<-"))) setGeneric("geneStudy<-",function(x,value) standardGeneric("geneStudy<-"))
if(is.null(getGeneric("coeff<-"))) setGeneric("coeff<-",function(x,value) standardGeneric("coeff<-"))
if(is.null(getGeneric("stdcoeff<-"))) setGeneric("stdcoeff<-",function(x,value) standardGeneric("stdcoeff<-"))
if(is.null(getGeneric("zscore<-"))) setGeneric("zscore<-",function(x,value) standardGeneric("zscore<-"))



### methods
### accessor functions
setMethod("exprs","mergeExprSet",function(object) return(object@data))
setMethod("geneNames","mergeExprSet",function(object) return(rownames(object@geneStudy)))
setMethod("names","mergeExprSet",function(x) return(names(x@data)))
setMethod("notes","mergeExprSet",function(object) return(object@notes))
setMethod("phenoData","mergeExprSet", function(object){
           x=exprs(object)
           y=list()
           for(i in 1:length(x)){ y[[i]]=pData(x[[i]])}
           return(y)})
setMethod("length","mergeExprSet",function(x) return(length(x@data)))
setMethod("geneStudy","mergeExprSet",function(x) return(x@geneStudy))
setMethod("coeff","mergeCoeff",function(x){
           y<-list()
           y[[1]]<-x@coeff
           y[[2]]<-"coeff"
           y[[3]]<-x@method
           return(y)})
setMethod("stdcoeff","mergeCoeff",function(x){
           y<-list()
           y[[1]]<-x@coeff.std
           y[[2]]<-"coeff.std"
           y[[3]]<-x@method
           return(y)})
setMethod("zscore","mergeCoeff",function(x){
           y<-list()
           y[[1]]<-x@zscore
           y[[2]]<-"zscore"
           y[[3]]<-x@method
           return(y)}) 
setMethod("modelOutcome","mergeExprSet",function(x,outcome,outcome2=NULL,method=c("linear","logistic","cox"),...) return(.model.outcome(x=x,outcome=outcome,outcome2=outcome2,method=method,...))) 
setMethod("intCor","mergeExprSet",function(x, method = c("pearson", "kendall", "spearman"),...) return(.intcor(x=x,method=method,...))) 
setMethod("plot","list",function(x,y,...) return(.plot.mergeCoeff(x=x,y=y,...)))
setMethod("plot","mergeCor",function(x,y,...) return(.plot.mergeCor(x=x,y=y,...)))
setMethod("hist","mergeCor",function(x,...) return(.hist.mergeCor(x=x,...)))
setMethod("intcorDens","mergeExprSet",function(x,method,...) return(.dens.mergeExprSet(x=x,method=method,...)))
setMethod("summary","mergeExprSet",function(object,...)  return(.summary.mergeExprSet(object=object,...)))

setMethod("intersection","mergeExprSet", function(x){
           nn<-length(x)
           tid<-geneNames(x)
           geneid  <- rep(1,nrow(x@geneStudy))
           for(i in 1:nn){
            geneid  <- geneid*(x@geneStudy[,i])
           }
           geneuid  <- geneNames(x)[geneid==1]
           cnote <- rep(0,nn)
           note  <-""
           k<-0
           for(i in 1:nn){
             mmatches<-match(geneuid,geneNames(exprs(x)[[i]]))
             exprs(exprs(x)[[i]])<-exprs(exprs(x)[[i]])[mmatches,]
             if(i==1) {
              ee<-exprs(exprs(x)[[i]])
              cnote[i]  <- ncol(exprs(exprs(x)[[i]]))
              k<-k+cnote[i]
              note<-paste(note,names(x)[i],": Column",1,"~Column",cnote[i],sep="")
             }
             else  {
              ee<-cbind(ee,exprs(exprs(x)[[i]]))
              cnote[i]  <- ncol(exprs(exprs(x)[[i]]))
              k<-k+cnote[i]
              note<-paste(note,", ",names(x)[i],": Column",k-cnote[i]+1,"~Column",k,sep="")
             }
           }
           y=new("exprSet",exprs=ee)
           geneNames(y)=geneuid
           notes(y)=note
           return(y)})
setMethod("cors","mergeCor",function(x) return(x@cors))
setMethod("pairwise.cors","mergeCor",function(x) return(x@pairwise.cors))
setMethod("integrative.cors","mergeCor",function(x){
           cc<-x@pairwise.cors
           UID  <- rownames(cc)
           avg.cc<-apply(cc,1,mean)
           names(avg.cc)<-UID
           return(avg.cc)
           })


subsetmES <- function(x, i,copy = TRUE){ 
  if(length(i)==1){return(exprs(x)[[i]])} else{
data=exprs(x)[i]
    notes=notes(x)
    gs=geneStudy(x)
    wh=(1:dim(gs)[1])[apply(gs[,i],1,sum)>0]
    gs=gs[wh,i]
return(new("mergeExprSet",data=data,geneStudy=gs,notes=notes))}}

setMethod("[", "mergeExprSet",
    function(x, i,j=NULL,drop=T) subsetmES(x, i, ...)
 )

### replacement functions
setReplaceMethod("exprs", "mergeExprSet", function(object, value){
 object@data<-value
object
})
setReplaceMethod("names", "mergeExprSet", function(x, value){
 names(x@data)<-value
x
})
setReplaceMethod("geneNames", "mergeExprSet", function(object, value){
nn<-length(object)
if(length(value)!=nrow(object@geneStudy)) stop("Replaced geneids should have the same length as the old genenames.") 
tid<-geneNames(object)
for(i in 1:nn){
 mmatches<-match(geneNames(exprs(object)[[i]]),tid)
 geneNames(exprs(object)[[i]])<-value[mmatches]
 
 idy  <- value[mmatches]
 y.avg  <- AverageDuplicates(exprs(exprs(object)[[i]]),idy)

 exprs(exprs(object)[[i]])  <- as.matrix(y.avg$data)
 geneNames(exprs(object)[[i]])  <- y.avg$acc
}
rownames(object@geneStudy)<-value
object
})
setReplaceMethod("notes", "mergeExprSet", function(object, value){
 object@notes<-value
object
})
setReplaceMethod("geneStudy", "mergeExprSet", function(x, value){
 x@geneStudy<-value
x})

setReplaceMethod("phenoData", "mergeExprSet", function(object, value){
nn<-length(object)
if (!is.element(class(value),"list")) stop("Replaced phenodata should be a list.") 
if(length(value)!=nn) stop("Replaced phenodata list should have the same length as the number of studies.") 
for(i in 1:nn){
 pData(exprs(object)[[i]])<-value[[i]]
}
object
})

setReplaceMethod("coeff", "mergeCoeff", function(x, value){
 x@coeff<-value
x})
setReplaceMethod("stdcoeff", "mergeCoeff", function(x, value){
 x@coeff.std<-value
x})
setReplaceMethod("zscore", "mergeCoeff", function(x, value){
 x@zscore<-value
x})

###subset mergeExprSet
check  <- function(x){
  if(!is.element(class(x),c("list","mergeExprSet","exprSet","matrix"))) stop("all data must be either a list, a mergeExprSet, a matrix or an exprSet")
}

mergeget  <- function(x){
  if (is.element(class(x),"mergeExprSet")){
   return(exprs(x))
  }
  if (is.element(class(x),"exprSet")){
   return(x)
  }
  if (is.element(class(x),"list")){
   if(length(x)!=4) stop("if you want to merge a list, a list should have at least four slots, 'expression matirx', 'phenodata', 'gene names' and 'notes'.")
   if (!is.matrix(x[[1]])) stop("first object of the input list must be an expression matrix.")
   if (!is.vector(x[[3]])) stop("third object of the input list must be a gene name vector.")
   if (is.null(x[[2]])) stop("second object of the input list can not be NULL.")
   tt  <- new("exprSet",exprs=x[[1]])
   pData(tt)  <- x[[2]]
   notes(tt)  <- x[[4]]
   geneNames(tt) <- x[[3]]
   return(tt)
  }
  if (is.element(class(x),"matrix")){
   if(is.null(rownames(x))) stop("if you want to merge matrix, rownames of matrix can not be NULL.")
   tt  <- new("exprSet",exprs=x)
   notes(tt)  <- ""
   geneNames(tt) <- rownames(x)
   return(tt)
  }
  stop("If you want to merge, the input object should be 'exprSet', 'list', 'matrix', or 'mergeExprSet'.")
}



AverageDuplicates  <- function(data.exprs,data.acc) {
  data.acc <- as.character(data.acc)
  dups<- rev(duplicated(rev(data.acc)))+duplicated(data.acc)
  dups<- ifelse(dups==0,0,1)
  
  if(sum(dups)>0) {
   data1.exprs<- data.exprs[dups==0,]
   data1.acc<-data.acc[dups==0]
   data2.exprs<- data.exprs[dups==1,]
   data2.acc<-data.acc[dups==1]
   data3.acc<-unique(data2.acc)
   data3.exprs<-matrix(NA,length(data3.acc),ncol(data.exprs))

   k <- 0
   for(i in data3.acc) {
    k <- k+1
    data3.exprs[k,]<-apply(data2.exprs[data2.acc==i,],2,mean,na.rm=TRUE)
   }
  
   data4.exprs<-rbind(data1.exprs,data3.exprs)
   data4.acc<-c(data1.acc,data3.acc)
   keep <- ifelse(as.character(data4.acc)=="",0,1)
 
   data.exprs <- data4.exprs[keep==1,]
   data.acc <- data4.acc[keep==1]
  }
  return(list(data=as.data.frame(data.exprs),acc=data.acc))
}


mergeExprs  <- function(...){
 arg  <- list(...)
 
 x  <- alist(...=)
 studynames<-list()
 k <- 0

 for(i in 1:length(arg)){
  check(arg[[i]])
  if (is.element(class(arg[[i]]),"mergeExprSet")){
   mm  <- mergeget(arg[[i]])
   studynames[[i]]<-names(arg[[i]])
   for(j in 1:length(arg[[i]])){
    k  <-k+1
    x[[k]]  <- mm[[j]]
   }
   
  }
  else {
   k  <-k+1
   x[[k]]  <- mergeget(arg[[i]])
   
   studynames[[i]]<-as.character(as.list(substitute(list(...)))[[i+1]])
  }
 }

 studynames<-unlist(studynames)
 tt  <- length(x)

 nnote  <- matrix(NA,tt,2)

 for (i in 1:tt){
   if(i==1) iid  <-as.matrix(geneNames(x[[i]]))
   else iid  <- rbind(iid,as.matrix(geneNames(x[[i]])))
   nnote[i,2]  <- notes(x[[i]])
 }

 iid  <- as.vector(sort(unique(iid)))
 
# generate the matrices with missing value "NA"

 for (i in 1:tt){
  y  <- as.matrix(exprs(x[[i]]))
  idy  <- geneNames(x[[i]])

  y.avg  <- AverageDuplicates(y,idy)

  exprs(x[[i]])  <- as.matrix(y.avg$data)
  geneNames(x[[i]])  <- y.avg$acc

  rownames(exprs(x[[i]]))  <- geneNames(x[[i]])
 }

# generate the vector with common id "1", o.w. "0"

 idmatrix  <- matrix(0,length(iid),tt)
 index  <- as.vector(nnote[,2])
 for (i in 1:tt){
  idx  <- geneNames(x[[i]])
  cc  <- match(iid,idx)
  idmatrix[,i]  <- ifelse(is.na(cc),0,1)
 }

 colnames(idmatrix)  <- studynames
 rownames(idmatrix)  <- iid
 names(x)  <- studynames

# generate the list that we want
 merged<-new("mergeExprSet",data=x,geneStudy=idmatrix,notes="")
 return(merged)
}





check.length  <- function(x,wh){
 if(length(pData(x)[,wh])!=ncol(exprs(x))) stop("Phenodata error: phenodata should have the same length as the number of columns of expression data.")
 return()
}

.model.outcome  <-function(x,method=NULL,outcome=NULL,outcome2=NULL){
 if(is.null(method)) stop("Specify the method you want to use.") 
 nn<-length(x)
 if (nn<=1) stop("Number of studies in the mergeExprSet should not less than 2.")
 
 if(method=="linear"|method=="logistic"){
  if(!is.null(outcome)){
   if(!is.element(class(outcome),"list")&!is.vector(outcome)) stop("Phenodata error: the data type of phenodata is not correct, please look at the help files.")
   if (is.element(class(outcome),"list")){
    if(length(outcome)!=nn)  stop("Phenodata error: the length of input phenodata list should be equal to the number of studies.")
    for(i in 1:nn){
     if(length(outcome[[i]])!=ncol(exprs(exprs(x)[[i]])))  stop("Phenodata error: phenodata should have the same length as the number of columns of expression data.")
    }
    out  <- outcome
   }
   if (is.vector(outcome)&!is.list(outcome)){
    if(length(outcome)!=nn) stop("Phenodata error: You should specify the phenodata.")
    else{
     out  <- alist(...=)
     for(i in 1:nn){
      check.length(exprs(x)[[i]],outcome[i])
      out[[i]]  <-  pData(exprs(x)[[i]])[,outcome[i]]
     }    
    }
   }
  }
  else{
    stop("Phenodata error: You should specify the phenodata.")
  }
 }

 if(method=="cox"){
  if(!is.null(outcome)&!is.null(outcome2)){
   if((!is.element(class(outcome),"list")&!is.vector(outcome))|(!is.element(class(outcome2),"list")&!is.vector(outcome2))) stop("Phenodata error: the data type of phenodata is not correct, please look at the help files.")
   if (is.element(class(outcome),"list")){
    if(length(outcome)!=nn)  stop("Phenodata error: the length of input phenodata list should be equal to the number of studies.")
    for(i in 1:nn){
     if(length(outcome[[i]])!=ncol(exprs(exprs(x)[[i]])))  stop("Phenodata error: phenodata should have the same length as the number of columns of expression data.")
    }
    out  <- outcome
   }

   if (is.vector(outcome)&!is.list(outcome)){
    if(length(outcome)!=nn) stop("Phenodata error: You should specify the phenodata.")
    else{
     out  <- alist(...=)
     for(i in 1:nn){
      check.length(exprs(x)[[i]],outcome[i])
      out[[i]]  <-  pData(exprs(x)[[i]])[,outcome[i]]
     }    
    }
   }
   
   if (is.element(class(outcome2),"list")){
    if(length(outcome2)!=nn)  stop("Phenodata error: the length of input phenodata list should be equal to the number of studies.")
    for(i in 1:nn){
     if(length(outcome2[[i]])!=ncol(exprs(exprs(x)[[i]])))  stop("Phenodata error: phenodata should have the same length as the number of columns of expression data.")
    }
    out2  <- outcome2
   }
   if (is.vector(outcome2)&!is.list(outcome)){
    if(length(outcome2)!=nn) stop("Phenodata error: You should specify the phenodata.")
    else{
     out2  <- alist(...=)
     for(i in 1:nn){
      check.length(exprs(x)[[i]],outcome2[i])
      out2[[i]]  <-  pData(exprs(x)[[i]])[,outcome2[i]]
     }    
    }
   }
  }
  else{
    stop("Phenodata error: You should specify the phenodata.")
  }
 }

 if(method!="linear"&method!="logistic"&method!="cox") stop("ARG should be one of linear, cox, logistic")
 
 geneid  <- rep(1,nrow(x@geneStudy))
 for(i in 1:nn){
  geneid  <- geneid*(x@geneStudy[,i])
 }
 geneuid  <- geneNames(x)[geneid==1]
 nuid  <- length(geneuid)

 beta  <- matrix(0,nuid,nn)
 stdbeta  <- matrix(0,nuid,nn)
 zscore  <- matrix(0,nuid,nn)

 for(i in 1:nn){
  matches1<-match(geneuid,geneNames(exprs(x)[[i]]))
  exprs1  <- exprs(exprs(x)[[i]])[matches1,]
  if(method=="linear"){
   outcome1  <- out[[i]]
   result1  <- apply(exprs1,1,.mergemodel,outcome=outcome1,method="linear")
   for(gene in 1:nuid){
    beta[gene,i]<-result1[[gene]][[1]]
    stdbeta[gene,i]<-result1[[gene]][[2]]
    zscore[gene,i]<-result1[[gene]][[3]]
   }
  }
  if(method=="logistic"){
   event1  <- out[[i]]
   result1  <- apply(exprs1,1,.mergemodel,event=event1,method="logistic")
   for(gene in 1:nuid){
    beta[gene,i]<-result1[[gene]][[1]]
    stdbeta[gene,i]<-result1[[gene]][[2]]
    zscore[gene,i]<-result1[[gene]][[3]]
   }
  }
  if(method=="cox"){
   outcome1  <- out[[i]]
   event1  <- out2[[i]]
   result1  <- apply(exprs1,1,.mergemodel,outcome=outcome1,event=event1,method="cox")
   for(gene in 1:nuid){
    beta[gene,i]<-result1[[gene]][[1]]
    stdbeta[gene,i]<-result1[[gene]][[2]]
    zscore[gene,i]<-result1[[gene]][[3]]
   }
  }
 }
 rownames(beta)<-rownames(stdbeta)<-rownames(zscore)<-geneuid
 
 nnote  <- rep(NA,nn)

 for(i in 1:nn){
  nnote[i]  <-notes(exprs(x)[[i]])
  if(is.null(notes(exprs(x)[[i]]))|is.na(notes(exprs(x)[[i]]))|notes(exprs(x)[[i]])=="") nnote[i]<-paste("study",i,sep=" ")
 }

 
 colnames(beta)<-colnames(stdbeta)<-colnames(zscore)<-names(x)

 result  <- new("mergeCoeff",coeff=beta, coeff.std=stdbeta,zscore=zscore,method=method)
 return(result)
}

.mergemodel <- function(x,outcome=NULL,event=NULL,method){
 if(method=="linear"){
  reg  <- lm(x~outcome,na.action=na.omit)
  sx  <- summary(reg)$sigma
  xx  <- x/sx
  stdbeta  <-  lm(outcome~xx,na.action=na.omit)$coeff[2]

  tmp  <-lm(outcome~x)
  beta  <- tmp$coeff[2]
  tvalue  <- summary(tmp)$coeff[2,3]

  return(list(beta=beta,stdbeta=stdbeta,tvalue=tvalue))
 }
 if(method=="cox"){
  ddtt  <- list(time = outcome,
                status = event,
                expr = x)
  expr  <- x
  status  <- event
  time  <- outcome

  cox.out  <- coxph( Surv(time,status) ~ expr, data=ddtt, na.action = na.omit)
  
  beta  <- cox.out$coeff
  
  tmp  <- status[!is.na(status)&!is.na(outcome)&!is.na(x)]
  nobs = sum(tmp==1)
  zscore  <- cox.out$coeff / sqrt(cox.out$var)

  std  <- list(time = outcome,
               status = event,
               expr = (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))
  cox.stdout  <- coxph( Surv(time,status) ~ expr, data=std, na.action = na.omit ) 
  stdbeta  <- cox.stdout$coeff
  
  return(list(beta=beta,stdbeta=stdbeta,zscore=zscore))
 }
 if(method=="logistic"){
  reg  <- glm(event~x,family=binomial,na.action = na.omit)
  beta  <- summary(reg)$coefficients[2,1]
  zscore  <- summary(reg)$coefficients[2,3]

  tmp  <- event[!is.na(x)&!is.na(event)]
  
  a  <- x[tmp==0]  
  b  <- x[tmp==1]
  na  <- sum(tmp==1)
  nb  <- length(tmp)-na
  
  if(na==0 || nb==0)  sp  <- sd(x)
  else  sp  <- (na*(sd(a,na.rm=TRUE)^2)+nb*(sd(b,na.rm=TRUE)^2))/length(x)

  stdx  <- x/sqrt(sp)
  regstd  <- glm(event~stdx,family=binomial,na.action = na.omit)

  stdbeta  <- summary(regstd)$coefficients[2,1]

  return(list(beta=beta,stdbeta=stdbeta,zscore=zscore)) }

}


.intcor <- function(x,method){
 nn  <- length(x)
 if (nn<=1) stop("Number of studies in the mergeExprSet should not less than 2.")
 
 geneid  <- rep(1,nrow(x@geneStudy))
 for(i in 1:nn){
  geneid  <- geneid*(x@geneStudy[,i])
 }
 geneuid  <- geneNames(x)[geneid==1]
 nuid  <- length(geneuid)

 pcor  <- alist(...=)
 nnote  <- rep(NA,nn)

 for(i in 1:nn){
  matches1<-match(geneuid,geneNames(exprs(x)[[i]]))
  exprs1  <- exprs(exprs(x)[[i]])[matches1,]
  
  pcor[[i]]  <- cor(t(exprs1),use="pairwise.complete.obs",method=method)
  nnote[i]  <-notes(exprs(x)[[i]])
  if(is.null(notes(exprs(x)[[i]]))|is.na(notes(exprs(x)[[i]]))|notes(exprs(x)[[i]])=="") nnote[i]<-paste("study",i,sep=" ")
 }
 
 names(pcor)  <- names(x)
 np<-nn*(nn-1)/2
 ppair  <- matrix(0,np,2)
 k  <- 1
 for(i in 1:(nn-1)){
  for(j in (i+1):nn){
   ppair[k,]<-c(i,j)
   k<-k+1
  }
 }
 
 icor<-matrix(NA,nuid,np)
 for(i in 1:np){
  CC1  <- pcor[[ppair[i,1]]]
  CC2  <- pcor[[ppair[i,2]]]

  for (gene in 1:nuid){
   c1  <- CC1[gene,][-gene]
   c2  <- CC2[gene,][-gene]
   icor[gene,i] <- cor(as.vector(c1),as.vector(c2),method=method)
  }
 }

 rownames(icor)<-geneuid
 result  <- new("mergeCor",cors=pcor, pairwise.cors=icor)
 return(result)

}
 


.plot.mergeCoeff  <-function(x,y,main=NULL,oma=NULL,...){
  if(is.null(main)){
   if(x[[2]]=="coeff") cc<-NULL
   if(x[[2]]=="coeff.std") cc<-"standardized coefficients"
   if(x[[2]]=="zscore") cc<-"zscore"
   if(x[[3]]=="linear") method<-"Linear Regression"
   if(x[[3]]=="cox") method<-"Cox Hazard Rate"
   if(x[[3]]=="logistic") method<-"Logistic Regression"
   main=paste(method,"Coefficients\n\n",cc,sep=" ")
   oma=c(2,2,8,2)
 }

 if(ncol(x[[1]])==2) plot(x[[1]][,1],x[[1]][,2],main=main,...)
 if(ncol(x[[1]])>2)  pairs(x[[1]],main=main,oma=oma,...)
 
 return()
}


.plot.mergeCor  <- function(x,y,labels=NULL,geneid=NULL,xlab=NA,ylab=NA,title=NULL,main=NULL,...){
 nn  <- length(x@cors)
 if(is.null(labels)) labels<-names(x@cors)
 if(length(labels)!=nn) stop("You need to specify the names for the studies.")
 if(nn<=1)  stop("You need to specify more than one study.")
 if(nn>2){
  np<-nn*(nn-1)/2 
  ppair  <- matrix(0,np,2)
  k  <- 1
  for(i in 1:(nn-1)){
   for(j in (i+1):nn){
    ppair[k,]<-c(i,j)
    k<-k+1
   }
  }
  indextmp<-ppair[,1]*nn+ppair[,2]
 
  ppair1  <- matrix(0,np,2)
  k  <- 1
  for(i in 1:(nn-1)){
   for(j in (i+1):nn){
    ppair1[k,]<-c(nn+1-i,nn+1-j)
    k<-k+1
   }
  }

  cc<-x@pairwise.cors
  cp<-x@cors
  k<-m<-1
  par(mfrow=c(nn,nn),oma=c(0,0,8,0))
  UID  <- rownames(cc)
  avg.cc<-apply(cc,1,mean)
  if(!is.null(geneid)){
   mmatches  <- match(geneid,UID)
   if(is.na(mmatches)) stop("No such common gene id.") 
  }
  else{
   mmatches<- (1:length(UID))[avg.cc==max(avg.cc)]
   score<-avg.cc[mmatches]
  }


  for(i in 1:nn){
   for(j in 1:nn){
    if(k<=np){
     if(ppair[k,1]==i&ppair[k,2]==j){
      cx<-cp[[ppair[k,1]]][mmatches,][-mmatches]
      cy<-cp[[ppair[k,2]]][mmatches,][-mmatches]
      plot(as.vector(cy),as.vector(cx),xlab=xlab,ylab=ylab,main=main,...)
      abline(h=0); abline(v=0)
      k=k+1
     }
    }
    if(m<=np){   
     if(ppair1[np+1-m,1]==i&ppair1[np+1-m,2]==j){
      pp1<-j*nn+i
      pp2<-c(1:np)[match(as.character(pp1),as.character(indextmp))]
      score  <- cc[,pp2][mmatches]
      aa<-c(.5,.5)
      bb<-c(.7,.8) 
      plot(aa,bb,xlim=c(0,1),ylim=c(0,1),xlab=NA,ylab=NA,pch=" ",xaxt="n",yaxt="n")
      text(.47,0.5,cex=1.2,paste("CORR = ",as.character(signif(score,digit=3))),col=4)
      m<-m+1
     }
    }
    if(i==j){
     aa<-c(.5,.5)
     bb<-c(.7,.8) 
     plot(aa,bb,xlim=c(0,1),ylim=c(0,1),xlab=NA,ylab=NA,pch=" ",xaxt="n",yaxt="n")
     text(.5,0.5,cex=1.2,paste(labels[i]),col=1)
    }
   }
  }
  score  <- avg.cc[mmatches]
  if(is.null(title)) title<-paste("Integrated Correlation\n\nGene",UID[mmatches],": average integrated score is",as.character(signif(score,digit=3)),sep=" ")
  mtext(title, line =0.5, cex=1.2,outer = TRUE)
  par(mfrow=c(1,1))
 }
 if(nn==2){
  cc<-x@pairwise.cors
  cp<-x@cors
  UID  <- rownames(cc)
  if(!is.null(geneid)){
   mmatches  <- match(geneid,UID)
   if(is.na(mmatches)) stop("No such common gene id.") 
  }
  if(is.null(geneid)) mmatches<- (1:length(UID))[cc==max(cc)]
  if(is.null(main)) main<-paste("Integrated Correlation\n\nGene",UID[mmatches])
  cx<-cp[[1]][mmatches,][-mmatches]
  cy<-cp[[2]][mmatches,][-mmatches]
  if(is.null(geneid)) score<-max(cc)
  else score<-cc[mmatches]
  par(oma=c(0,0,0,3))
  if(is.na(xlab)) plot(as.vector(cx),as.vector(cy),ylab=ylab,xlab=paste("CORR = ",as.character(signif(score,digit=3))),main=main,...)
  else plot(as.vector(cx),as.vector(cy),ylab=ylab,xlab=paste("\n",xlab,"\nCORR = ",as.character(signif(score,digit=3))),main=main,...)
  abline(h=0); abline(v=0)
 }
 return()
}


.hist.mergeCor  <- function(x,labels=NULL,main=NULL,xlab=NULL,title=NULL,...){
 nn  <- length(x@cors)
 if(is.null(labels)) labels<-names(x@cors)
 if(length(labels)!=nn) stop("You need to specify the names for the studies.")
 if(nn<=1)  stop("You need to specify more than one study.")

 if(is.null(main)) main<-NA
 if(is.null(xlab)) xlab<-NA
 if(nn>2){
  np<-nn*(nn-1)/2 
  ppair  <- matrix(0,np,2)
  k  <- 1
  for(i in 1:(nn-1)){
   for(j in (i+1):nn){
    ppair[k,]<-c(i,j)
    k<-k+1
   }
  }
  indextmp<-ppair[,1]*nn+ppair[,2]  

  ppair1  <- matrix(0,np,2)
  k  <- 1
  for(i in 1:(nn-1)){
   for(j in (i+1):nn){
    ppair1[k,]<-c(nn+1-i,nn+1-j)
    k<-k+1
   }
  }

  cp<-x@pairwise.cors
  k<-m<-1
  par(mfrow=c(nn,nn),oma=c(0,0,5,0))
  
  
  for(i in 1:nn){
   for(j in 1:nn){
    if(k<=np){
     if(ppair[k,1]==i&ppair[k,2]==j){
      hist(cp[,k],xlab=xlab,main=main,...)
      k=k+1
     }
    }
    if(m<=np){   
     if(ppair1[np+1-m,1]==i&ppair1[np+1-m,2]==j){
      pp1<-j*nn+i
      pp2<-c(1:np)[match(as.character(pp1),as.character(indextmp))]
      hist(cp[,pp2],xlab=xlab,main=main,...)
      m<-m+1
     }
    }
    if(i==j){
     aa<-c(.5,.5)
     bb<-c(.7,.8) 
     plot(aa,bb,xlim=c(0,1),ylim=c(0,1),xlab=NA,ylab=NA,pch=" ",xaxt="n",yaxt="n")
     text(.5,0.5,cex=1.2,paste(labels[i]),col=1)
    }
   }
  }
  if(is.null(title)) title<-"Integrated Correlation"
  mtext(title, line =0.5, cex=1.2,outer = TRUE)
  par(mfrow=c(1,1))
 }
 if(nn==2){
  if(is.na(main)) main<-"Integrated Correlation"
  cp<-x@pairwise.cors
  hist(cp,main=main,xlab=xlab,...)
 }

 return()
}


.dens.mergeExprSet  <- function(x,main=NA,x.legend=NULL,y.legend=NULL,cex.legend=NULL,title=NULL,method=NULL,...){
 nn  <- length(x)
 if (nn<=1) stop("Number of studies in the mergeExprSet should not less than 2.")

 nnote  <- rep(NA,nn)

 for(i in 1:nn){
  nnote[i]  <-notes(exprs(x)[[i]])
  if(is.null(notes(exprs(x)[[i]]))|is.na(notes(exprs(x)[[i]]))|notes(exprs(x)[[i]])=="") nnote[i]<-paste("study",i,sep=" ")
 }

 geneid  <- rep(1,nrow(x@geneStudy))
 for(i in 1:nn){
  geneid  <- geneid*(x@geneStudy[,i])
 }
 geneuid  <- geneNames(x)[geneid==1]
 nuid  <- length(geneuid)

 if(nn>2){
  np<-nn*(nn-1)/2 
  ppair  <- matrix(0,np,2)
  k  <- 1
  for(i in 1:(nn-1)){
   for(j in (i+1):nn){
    ppair[k,]<-c(i,j)
    k<-k+1
   }
  }
  
  ppair1  <- matrix(0,np,2)
  k  <- 1
  for(i in 1:(nn-1)){
   for(j in (i+1):nn){
    ppair1[k,]<-c(nn+1-i,nn+1-j)
    k<-k+1
   }
  }


  exp1<-matrix(NA,nuid,np)
  exp2<-matrix(NA,nuid,np)
  obs<-matrix(NA,nuid,np)
  
  for(i in 1:np){
   matches1  <- match(geneuid,geneNames(exprs(x)[[ppair[i,1]]]))
   matches2  <- match(geneuid,geneNames(exprs(x)[[ppair[i,2]]]))
   exprs1  <- exprs(exprs(x)[[ppair[i,1]]])[matches1,]
   exprs2  <- exprs(exprs(x)[[ppair[i,2]]])[matches2,]
   exprs1[is.na(exprs1)]  <-  0
   exprs2[is.na(exprs2)]  <-  0

   cor1  <- cor(t(exprs1),use="pairwise.complete.obs",method=method)
   cor2  <- cor(t(exprs2),use="pairwise.complete.obs",method=method)
  
   for (gene in 1:nuid){
    c1  <- cor1[gene,][-gene]
    c2  <- cor2[gene,][-gene]
    obs[gene,i] <- cor(as.vector(c1),as.vector(c2),method=method)
   }

   exprstmp<-matrix(NA,nrow(exprs1),ncol(exprs1))
   

   for (gene in 1:nuid){
    subset <- sample(1:ncol(exprs1))
    exprstmp[gene,] <- exprs1[gene,subset]
   }

   cor1  <- cor(t(exprstmp),use="pairwise.complete.obs",method=method)
   
   for (gene in 1:nuid){
    c1  <- cor1[gene,][-gene]
    c2  <- cor2[gene,][-gene]
    exp1[gene,i] <- cor(as.vector(c1),as.vector(c2),method=method)
   }

   for (gene in 1:nuid){
    subset <- sample(1:ncol(exprs1))
    exprstmp[gene,]  <- exprs1[gene,subset]
   }

   cor1  <- cor(t(exprstmp),use="pairwise.complete.obs",method=method)
   
   for (gene in 1:nuid){
    c1  <- cor1[gene,][-gene]
    c2  <- cor2[gene,][-gene]
    exp2[gene,i] <- cor(as.vector(c1),as.vector(c2),method=method)
   }
  }
     
  k<-m<-1
  par(mfrow=c(nn,nn),oma=c(0,0,8,0))
  
  for(i in 1:nn){
   for(j in 1:nn){
    if(k<=np){
     if(ppair[k,1]==i&ppair[k,2]==j){
      d1<-exp1[,k]
      d2<-exp2[,k]
      d<-obs[,k]
      xmax<-max(c(density(d1)$x,density(d2)$x,density(d)$x))
      xmin<-min(c(density(d1)$x,density(d2)$x,density(d)$x))
      ymax<-max(c(density(d1)$y,density(d2)$y,density(d)$y))
      range<-xmax-xmin
      plot(density(d),xlim=c(xmin-range/10,xmax+range/10),ylim=c(0,ymax+ymax/5),main=main,...)
      lines(density(d1),col=2,...)
      lines(density(d2),col=3,...)
     
      k=k+1
     }
    }
    if(m<=np){
     if(ppair1[np+1-m,1]==i&ppair1[np+1-m,2]==j){
      if(ppair1[np+1-m,1]==nn&ppair1[np+1-m,2]==1){
       aa<-c(.5,.5) 
       bb<-c(.7,.8) 
       plot(aa,bb,xlim=c(0,1),ylim=c(0,1),xlab=NA,ylab=NA,pch=" ",xaxt="n",yaxt="n")
       legend.text  <- c("observed","resampling 1","resampling 2")
       op <- par(bg="white")
       if(is.null(x.legend)) x.legend<-.1
       if(is.null(y.legend)) y.legend<-.9
       if(is.null(cex.legend)) cex.legend<-1
       legend(x.legend,y.legend, paste(legend.text), col=1:3 ,lty=1,cex=cex.legend)
      }
      else {
       aa<-c(.5,.5) 
       bb<-c(.7,.8) 
       plot(aa,bb,xlim=c(0,1),ylim=c(0,1),xlab=NA,ylab=NA,pch=" ",xaxt="n",yaxt="n")
      }
      m=m+1
     }
    }
    if(i==j){
     aa<-c(.5,.5) 
     bb<-c(.7,.8) 
     plot(aa,bb,xlim=c(0,1),ylim=c(0,1),xlab=NA,ylab=NA,pch=" ",xaxt="n",yaxt="n")
     text(.5,0.5,cex=1.2,paste(names(x)[i]),col=1)
    }
   }
  }
  if(is.null(title)) title<-"Integrated Correlation:\n\ntrue and null distributions"
  mtext(title, line =0.5, cex=1.2,  outer = TRUE)
  par(mfrow=c(1,1)) 
 }

 if(nn==2){
  geneid  <- x@geneStudy[,1]*x@geneStudy[,2]
 
  geneuid  <- geneNames(x)[geneid==1]
  nuid  <- length(geneuid)
  
  matches1  <- match(geneuid,geneNames(exprs(x)[[1]]))
  matches2  <- match(geneuid,geneNames(exprs(x)[[2]]))
  exprs1  <- exprs(exprs(x)[[1]])[matches1,]
  exprs2  <- exprs(exprs(x)[[2]])[matches2,]
  exprs1[is.na(exprs1)]  <-  0
  exprs2[is.na(exprs2)]  <-  0

  cor1  <- cor(t(exprs1),use="pairwise.complete.obs",method=method)
  cor2  <- cor(t(exprs2),use="pairwise.complete.obs",method=method)
  
  d<-rep(NA,nuid)
  d1<-rep(NA,nuid)
  d2<-rep(NA,nuid)
  for (gene in 1:nuid){
   c1  <- cor1[gene,][-gene]
   c2  <- cor2[gene,][-gene]
   d[gene] <- cor(as.vector(c1),as.vector(c2),method=method)
  }
  
  exprstmp<-matrix(NA,nrow(exprs1),ncol(exprs1))
      
  for (gene in 1:nuid){
   subset <- sample(1:ncol(exprs1))
   exprstmp[gene,] <- exprs1[gene,subset]
  }

  cor1  <- cor(t(exprstmp),use="pairwise.complete.obs",method=method)
   
  for (gene in 1:nuid){
   c1  <- cor1[gene,][-gene]
   c2  <- cor2[gene,][-gene]
   d1[gene] <- cor(as.vector(c1),as.vector(c2),method=method)
  }

  for (gene in 1:nuid){
   subset <- sample(1:ncol(exprs1))
   exprstmp[gene,] <- exprs1[gene,subset]
  }

  cor1  <- cor(t(exprstmp),use="pairwise.complete.obs",method=method)
   
  for (gene in 1:nuid){
   c1  <- cor1[gene,][-gene]
   c2  <- cor2[gene,][-gene]
   d2[gene] <- cor(as.vector(c1),as.vector(c2))
  } 
  
  xmax<-max(c(density(d1)$x,density(d2)$x,density(d)$x))
  xmin<-min(c(density(d1)$x,density(d2)$x,density(d)$x))
  ymax<-max(c(density(d1)$y,density(d2)$y,density(d)$y))
  range<-xmax-xmin
  if(is.null(main)) main<-"Integrated Correlation:\ntrue and null distributions"
  plot(density(d),xlim=c(xmin-range/10,xmax+range/10),ylim=c(0,ymax+ymax/5),main=main,...)
  lines(density(d1),col=2,...)
  lines(density(d2),col=3,...)   
  legend.text  <- c("observed","resampling 1","resampling 2")
  op <- par(bg="white")
  if(is.null(x.legend)) x.legend<-xmax-range/3
  if(is.null(y.legend)) y.legend<-ymax
  if(is.null(cex.legend)) cex.legend<-.9
  legend(x.legend,y.legend, paste(legend.text), col=1:3,lty=1,cex=cex.legend)
  
 }
 return()

}


.summary.mergeExprSet<-function(object){
 report  <- alist(...=)
 nn<-length(object)

 nnote  <- rep(NA,nn)

 for(i in 1:nn){
  nnote[i]  <-notes(exprs(object)[[i]])
  if(is.null(notes(exprs(object)[[i]]))|is.na(notes(exprs(object)[[i]]))|notes(exprs(object)[[i]])=="") nnote[i]<-paste("study",i,sep=" ")
 } 
 report[[1]]<-matrix(NA,2,nn)
 report[[1]][1,]<-names(object)
 for(i in 1:nn){
  report[[1]][2,i]<-nrow(exprs(exprs(object)[[i]]))
 }

 report[[2]]<-matrix(NA,2,nn)
 report[[2]][1,]<-names(object)
 for(i in 1:nn){
  report[[2]][2,i]<-ncol(exprs(exprs(object)[[i]]))
 }

 report[[3]]<-notes(object) 

 names(report)<-c("Number of Genes in Each Study","Number of Samples in Each Study","Notes")
 return(report)
}
