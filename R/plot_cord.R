#' Plots an ordination plot of latent variables and their corresponding coefficients (biplot).
#'
#' @param obj is a cord object, e.g. from output of \code{cord}
#' @param biplot \code{TRUE} if both latent variables and their coefficients are plotted, \code{FALSE} if only latent variables
#' @param site.col site number colour (default is black), vector of length equal to the number of sites 
#' @param sp.col species name colour (default is blue), vector of length equal to the number of sites (if arrow=TRUE)
#' @param alpha scaling factor for ratio of scores to loadings (default is 0.9)
#' @param arrow should arrows be plotted for species loadings (default is TRUE)
#' @param ...	other parameters to be passed through to plotting functions. 
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
#' plot(spid_lv,biplot = TRUE,site.col=cols)

plot.cord <- function(obj, biplot = FALSE,site.col="black",sp.col="blue",alpha=0.7,arrow=TRUE, ...) {
    
    if(dim(spid_lv$loadings)[2]==1){
      stop("this function does not plot a single factor")
    }
    if(dim(spid_lv$loadings)[2]!=2){
      warning("plotting first 2 factors")
    }
    #species labels from the original data
    labs = colnames(obj$obj$data$abund)
    
    #extract scores and loadings
    loadings=obj$loadings[,1:2]
    scores=t(obj$scores[1:2,])
  
    #calculate scaling factor
    alpha_plot=sqrt(max(apply(scores^2,1,max)))/sqrt(max(apply(loadings^2,1,max)))*alpha
    
    #plot graph
    plot(scores,pch=16,ylab="Latent variable 2",
         xlab="Latent variable 1", type='n',...)
    abline(h=0,col="gray")
    abline(v=0,col="gray")
    text(scores,label=1:nrow(scores),col=site.col)
    if(biplot){
      loadings=loadings*alpha_plot
      if(arrow){
        arrows(0,0,loadings[,1],loadings[,2],length=0,col="dark gray")
      }
      text(loadings[,1],loadings[,2],labels = labs,col=sp.col,cex = 0.8)
    }
    
}

