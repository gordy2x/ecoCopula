

#' Plot graph of species intearactions.
#'
#' @param obj is a cord object, e.g. from output of \code{cord}
#' @return an ordination plot.
#' @export
#' @examples
#' data(spider)
#' abund <- mvabund(spider$abund)
#' X <- spider$x
#' spider_mod=manyglm(abund~1)
#' spid_lv=cord(spider_mod)
#' #colour sites accoring to second column of x (bare sand)
#' cols=ifelse(spider$x[,2]>0,"black","red")
#' ordiplot(spid_lv,biplot = TRUE,col=cols)

ordiplot.cord = function(obj,biplot = FALSE, ...) {
    
    #species labels from the original data
    labs = colnames(obj$obj$data$abund)
    
    #extract scores and loadings
    loadings=obj$loadings[[1]][,1:2]
    scores=t(obj$scores[[1]][1:2,])
  
    #calculate scaling factor
    alpha=sqrt(max(apply(scores^2,1,sum)))/sqrt(max(apply(loadings^2,1,sum)))*0.9
    
    #plot graph
    plot(scores,pch=16,ylab="Latent variable 2",
         xlab="Latent variable 1", type='n')
    abline(h=0,col="gray")
    abline(v=0,col="gray")
    text(scores,label=1:nrow(scores),...)
    if(biplot){
      loadings=loadings*alpha
      # arrows(0,0,loadings[,1],loadings[,2],length=0,col="dark gray")
      text(loadings[,1],loadings[,2],labels = labs,col="blue",cex = 0.7)
    }
    
}

