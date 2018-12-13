
#' Fit gaussian copula graphical lasso, chose lambda with BIC.
#'
#' @param obj object of class manyglm, or manyany with "clm"
#' @param lambda values of shrinkage parameter lambda for model selection (optional)
#' @param n.lambda number of values of shrinkage parameter lambda for model selection, ignored if lambda suplied
#' @param n.samp number of residuals used for importance sampling (optional)
#' @return best.graph a list with the paramaters for the 'best' graphical model, chosen by BIC
#' @return all.graphs a list with likelihood and BIC for all models along lambda path
#' @export 
#' @examples
#' data(spider)
#' abund <- mvabund(spider$abund)
#' X <- spider$x
#' spider.mod=manyglm(abund~X)
#' spid.graph=saint(spider.mod)
#' plot(spid.graph)
#' 
#' 
saint<-function(obj,lambda=NULL,n.lambda=100,n.samp=500,seed=1)
{
  
  if(!is.numeric(seed))
    stop("seed must be numeric")
  
  if(floor(n.samp)!=ceiling(n.samp))
    stop("n.samp must be an integer")
  
  if(floor(n.lambda)!=ceiling(n.lambda))
    stop("n.lambda must be an integer")
  
  if(!(is.numeric(lambda)|is.null(lambda)))
    stop("lambda must be numeric")
  
  if(!is.null(lambda))
    warning(" 'best' model selected among supplied lambda only")
  
  if(any(lambda<0))
    stop("lambda must be non negative")
  
  if(!class(obj)[1]%in% c("manyany","manyglm"))
    stop("please supply an manyglm, or manyany (clm) object")
  
  if(class(obj)[1]=="manyany" & class(obj)[2]!="clm")
    stop("saint function only supposts manyany with clm")
  
  #always same result unless specified otherwise
  set.seed(seed)
  
  #simulate full set of residulas n.samp times
  res=simulate.res.S(obj,n.res=n.samp)
  
  if(is.null(lambda)){
    #starting values for log10 lambda
    current=c(-16,seq(-10,10,length.out=20))

    #proportion of non-zero cond indep
    k.frac=full.graph.many(obj,10^current,res)$k.frac
    #which lambdas give full and empty matrices
    new.min=max(which(k.frac==1))
    new.max=min(which(k.frac==0))
    
    #now sequence of lambda values, between the two
    lambda=10^c(-16,seq(current[new.min],current[new.max],length.out=n.lambda))
  } else{
    n.lambda=length(lambda)
  }
  ag=full.graph.many(obj,lambda,res)
  
  k.frac=ag$k.frac
  BIC.graph=ag$BIC
  logL=ag$logL
  
  
  best=min(which(BIC.graph==min(BIC.graph)))
  Th.best=ag$Th.out[[best]]
  Sig.best=ag$Sig.out[[best]]
  
  graph.out=as.matrix((Th.best!=0)*1)
  best.graph=list(graph=graph.out,prec=Th.best,cov=Sig.best,Y=obj$y,logL=logL[[best]],sparsity=k.frac[best])
  all.graphs=list(lambda.opt=lambda[best],logL=logL,BIC=BIC.graph,lambda=lambda,
                  k.frac=k.frac)
  out=list(best.graph=best.graph,all.graphs=all.graphs,obj=obj)
  class(out)="saint"
  return(out)
  
}