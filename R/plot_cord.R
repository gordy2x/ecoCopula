#' Plots an ordination of latent variables and their corresponding coefficients (biplot).
#'
#' @param x is a cord object, e.g. from output of \code{cord}
#' @param biplot \code{TRUE} if both latent variables and their coefficients are plotted, \code{FALSE} if only latent variables
#' @param site.col site number colour (default is black), vector of length equal to the number of sites 
#' @param sp.col species name colour (default is blue), vector of length equal to the number of sites (if arrow=TRUE)
#' @param alpha scaling factor for ratio of scores to loadings (default is 0.7)
#' @param arrow should arrows be plotted for species loadings (default is TRUE)
#' @param site.text should sites be labeled by rown names of data (default is FALSE, points are drawn)
#' @param labels the labels for sites and species (for biplots only) (default is data labels)
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
#' plot(spid_lv,biplot = FALSE,site.col=cols, site.text = TRUE)

plot.cord <- function(x, biplot = FALSE,site.col="black",sp.col="blue",
                      alpha=0.7,arrow=TRUE, site.text=FALSE,
                      labels=dimnames(x$obj$fitted),...) {
    
    if(dim(x$loadings)[2]==1){
      stop("this function does not plot a single factor")
    }
    if(dim(x$loadings)[2]!=2){
      warning("plotting first 2 factors")
    }
  
    #extract scores and loadings
    loadings=x$loadings[,1:2]
    scores=x$scores[,1:2]
  
    #calculate scaling factor
    alpha_plot=sqrt(max(apply(scores^2,1,max)))/sqrt(max(apply(loadings^2,1,max)))*alpha
    
    #plot graph
    plot(scores,ylab="Latent variable 2",
         xlab="Latent variable 1", type='n',...)
    abline(h=0,col="gray")
    abline(v=0,col="gray")
    
    if(site.text){
      text(scores,label=labels[[1]],col=site.col)
    }else{
      points(scores,col=site.col,pch=16)
    }
    if(biplot){
      loadings=loadings*alpha_plot
      if(arrow){
        arrows(0,0,loadings[,1],loadings[,2],length=0,col="dark gray")
      }
      text(loadings[,1],loadings[,2],labels = labels[[2]],col=sp.col,cex = 0.8)
    }
    
}

