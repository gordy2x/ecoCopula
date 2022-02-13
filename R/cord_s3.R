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
#' @importFrom grDevices devAskNewPage
#' @export
#' @examples
#' X <- spider$x
#' abund <- spider$abund
#' spider_mod <- stackedsdm(abund,~1, data = X, ncores=2) 
#' spid_lv=cord(spider_mod)
#' #colour sites according to second column of x (bare sand)
#' cols=ifelse(spider$x[,2]>0,"black","red")
#' plot(spid_lv,biplot = FALSE,site.col=cols, site.text = TRUE)
#' 
#'\donttest{
#'library(ggplot2)
#'library(RColorBrewer)
#'alpha= 2.5
#'site_res <- data.frame(spid_lv$scores,X)
#'sp_res <- data.frame(spid_lv$loadings,species=colnames(abund))
#'ggplot()+
#'  geom_point(aes(x=Factor1,y=Factor2,color = reflection ),site_res)+
#'  geom_text(aes(x = Factor1*alpha, y = Factor2*alpha,label = species),data=sp_res)+
#'  scale_color_gradientn(colours = brewer.pal(n = 10, name = "PuOr"))+
#'  theme_classic()
#'}
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
    
    #don't ask for new plot but reset globally after
    oask <- devAskNewPage()
    devAskNewPage(FALSE)
    on.exit(devAskNewPage(oask))
    
    #plot graph
    plot(scores,ylab="Latent variable 2",
         xlab="Latent variable 1", type='n',...)
    abline(h=0,col="gray")
    abline(v=0,col="gray")
    
    if(site.text){
      text(scores,label=labels[[1]],col=site.col)
    }else{
      graphics::points(scores,col=site.col,pch=16)
    }
    if(biplot){
      loadings=loadings*alpha_plot
      if(arrow){
        arrows(0,0,loadings[,1],loadings[,2],length=0,col="dark gray")
      }
      text(loadings[,1],loadings[,2],labels = labels[[2]],col=sp.col,cex = 0.8)
    }
    
}



#' Print function for cord object
#'
#' @param x is a cord object, e.g. from output of \code{\link{cord}}.
#' @param ... not used
#' @seealso \code{\link{cord}}
#' @export
#' @examples
#' abund <- spider$abund
#' spider_mod <- stackedsdm(abund,~1, data = spider$x, ncores=2) 
#' spid_lv=cord(spider_mod)
#' print(spid_lv)
print.cord = function (x, ...) 
{
  cat("\nCall:\n", paste(deparse(x$obj$call), sep = "\n", 
                         collapse = "\n"), "\n\n", sep = "")
  
  cat("Pairwise associations:\n")
  print(paste0(ncol(x$loadings)," latent variables"))
  
  cat("\n")
  invisible(x)
}

#' Summary function for cgr object
#'
#' @param object is a cord object, e.g. from output of \code{\link{cgr}}.
#' @param ... not used
#' @seealso \code{\link{cord}}
#' @export
#' @examples
#' abund <- spider$abund[,1:5]
#' spider_mod <- stackedsdm(abund,~1, data = spider$x, ncores=2) 
#' spid_lv=cord(spider_mod)
#' summary(spid_lv)
summary.cord = function (object, ...) 
{
  cat("\nCall:\n", paste(deparse(object$obj$call), sep = "\n", 
                         collapse = "\n"), "\n\n", sep = "")
  
  cat("Loadings:\n")
  print(round(object$loadings,3))
  
  cat("\n")
  invisible(object)
}


#'Simulates new data from a given cord object
#'
#' @param object is a cord object, e.g. from output of \code{cord}
#' @param nsim Number of simulations, defaults to 1. If nsim > 1, the simulated data will be
#' appended.
#' @param seed Random number seed, defaults to a random seed number.
#' @param newdata A data frame in which to look for X covariates with which to simulate.
#' @param ... not used
#' Defaults to the X covariates in the fitted model.
#' @examples
#' abund = spider$abund
#'
#' spider_mod_ssdm = stackedsdm(abund,~1, data = spider$x, ncores=2)
#' spid_lv_ssdm = cord(spider_mod_ssdm)
#' simulate(spid_lv_ssdm, nsim=2)
#' 
#'\donttest{
#' # using mvabund
#' library(mvabund) #for manyglm
#' abund=mvabund(abund)
#' spider_mod = manyglm(abund~1)
#' spid_lv = cord(spider_mod)
#' simulate(spid_lv)
#'
#' spider_mod_X = manyglm(abund ~ soil.dry + bare.sand, data=spider$x)
#' spid_lv_X = cord(spider_mod_X)
#' Xnew = spider$x[1:10,]
#' simulate(spid_lv_X, newdata = Xnew)
#' simulate(spid_lv_X, nsim=2, newdata = Xnew)
#'
#' spider_mod_X_ssdm = stackedsdm(abund, formula_X = ~. -bare.sand, data = spider$x, ncores=2)
#' spid_lv_X_ssdm = cord(spider_mod_X_ssdm)
#' simulate(spid_lv_X_ssdm, newdata = Xnew)
#' }
#' @importFrom stats simulate
#' @export
simulate.cord = function(object, nsim=1, seed=NULL, newdata=object$obj$data, ...) {
  
  # code chunk from simulate.lm to select seed
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  if (is.null(seed)) { 
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  check_object_family(object)
  newdata = reshape_newdata(object, nsim, newdata)
  prs = suppressWarnings(
    predict(object$obj, type = "response", newdata = newdata)
  ) # warning for family=poisson suppressed
  newY = simulate_newY(object, prs)
  
  return (newY)
}

#' Simulate new Y
#'
#' Simulates new Y
#'
#' @param object 
#' @param prs 
#'
#' @noRd
simulate_newY = function(object, prs) {
  nRow = nrow(prs)
  nVar = ncol(prs)
  newY = matrix(NA, nRow, nVar)
  
  # simulate multivariate normal random variables
  sim = MASS::mvrnorm(nRow, mu = rep(0, times = nVar), object$sigma)
  
  # turn simulated variables into abundances
  if (as.character(object$obj$call[[1]]) %in% c("manyglm", "stackedsdm")) {
    for (iVar in 1:nVar) {
      if (all(object$obj$family == "negative.binomial")) {
        size = get_size(object)
        newY[,iVar] = qnbinom(pnorm(sim[,iVar]), size = size[iVar], mu = prs[,iVar])
      } else if (object$obj$call$family == "poisson") {
        newY[,iVar] = qpois(pnorm(sim[,iVar]), lambda = prs[,iVar])
      } else if (object$obj$call$family == "binomial") {
        newY[,iVar] = qbinom(pnorm(sim[,iVar]), size = 1, prob = prs[,iVar])
      } else {
        stop("'family'", object$obj$family, "not recognized")
      }
    }
  } else if (object$obj$call[[1]] == "manylm") {
    df.residual = object$obj$df.residual
    sigma2 = apply((object$obj$y - object$obj$fitted)^2, 2, sum)/df.residual
    for (iVar in 1:nVar) {
      newY[,iVar] = sim[,iVar] * sqrt(sigma2[iVar])
    }
    newY = newY + prs
  } else {
    stop("'class'", class(object$obj), "not supported")
  }
  
  colnames(newY) = colnames(object$obj$y)
  return (newY)
}
