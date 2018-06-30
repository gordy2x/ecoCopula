
#' Fit gaussian copula graphical lasso, chose lambda with BIC.
#'
#' @param manyglm.obj object of class manyglm.
#' @param n.lambda number of values of shrinkage parameter lambda for model selection.
#' @param n.samp number of residuals used for importance sampling.
#' @return best.graph a list with the paramaters for the best graphical model
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
saint<-function(manyglm.obj,n.lambda=100,n.samp=200)
{
  set.seed(1)
  current=c(-16,seq(-10,10,length.out=20))
  res=simulate.res.S(manyglm.obj,n.res=n.samp)
  k.frac=full.graph.many(manyglm.obj,10^current,res)$k.frac
  new.min=max(which(k.frac==1))
  new.max=min(which(k.frac==0))

  
  final=c(-16,seq(current[new.min],current[new.max],length.out=n.lambda))
  ag=full.graph.many(manyglm.obj,10^final,res)
  k.frac=ag$k.frac
  BIC.graph=ag$BIC
  logL=ag$logL
  
  
  best=min(which(BIC.graph==min(BIC.graph)))
  Th.best=ag$Th.out[[best]]
  Sig.best=ag$Sig.out[[best]]
  
  graph.out=as.matrix((Th.best!=0)*1)
  best.graph=list(graph=graph.out,prec=Th.best,cov=Sig.best,Y=manyglm.obj$y,logL=logL[[best]],sparsity=k.frac[best])
  all.graphs=list(lambda.opt=10^final[best],logL=logL,BIC=BIC.graph,lambda=10^final,
                  k.frac=k.frac)
  out=list(best.graph=best.graph,all.graphs=all.graphs,manyglm.obj=manyglm.obj)
  class(out)="saint"
  return(out)
  
}