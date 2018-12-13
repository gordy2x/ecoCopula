
#' Model based ordination with gaussian copulas
#'
#' @param obj object of class manyglm, or manyany with "clm"
#' @param nlv number of latent variables (default is 2)
#' @param n.samp number of residuals used for importance sampling (default is 500)
#' @return 
#' @export 
#' @examples
#' data(spider)
#' abund <- mvabund(spider$abund)
#' X <- spider$x
#' spider.mod=manyglm(abund~X)
#' spid.lv=cord(spider.mod)
#' plot(spid.lv)
#' 
#' 
cord<-function(obj,nlv=2,n.samp=500,seed=1)
{
  
  if(!is.numeric(seed))
    stop("seed must be numeric")
  
  if(floor(n.samp)!=ceiling(n.samp))
    stop("n.samp must be an integer")
  
  if(floor(nlv)!=ceiling(nlv))
    stop("nlv must be an integer")
  
  
  #always same result unless specified otherwise
  set.seed(seed)
 
  #simulate full set of residulas n.samp times
  res=simulate.res.S(obj,n.res=n.samp)
  
  #carry out factor analysis
  out=fact.many(obj,nlv,res)

  class(out)="cord"
  return(out)
  
}