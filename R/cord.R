#' Model based ordination with Gaussian copulas
#'
#' @param  obj object of either class \code{\link[mvabund]{manyglm}}, 
#' or  \code{\link[mvabund]{manyany}} with ordinal models \code{\link[ordinal]{clm}}
#' @param nlv number of latent variables (default = 2, for plotting on a scatterplot)
#' @param n.samp integer (default = 500), number of sets residuals used for importance sampling 
#' (optional, see detail)
#' @param seed integer (default = 1), seed for random number generation (optional, see detail)

#' @section Details:
#' \code{cord} is used to fit a Gaussian copula factor analytic model to multivariate discrete data, such as co-occurence (multi species) data in ecology. The model is estimated using importance sampling with \code{n.samp} sets of randomised quantile or "Dunn-Smyth" residuals (Dunn & Smyth 1996), and the \code{\link{factanal}} function. The seed is controlled so that models with the same data and different predictors can be compared.  
#' @return
#' \code{loadings} latent factor loadings
#' \code{scores} latent factor scores
#' \code{sigma} covariance matrix estimated with \code{nlv} latent variables
#' \code{theta} precision matrix estimated with \code{nlv} latent variables
#' \code{BIC} BIC of estimated model
#' \code{logL} log-likelihood of estimated model
#' @section Author(s):
#' Gordana Popovic <g.popovic@unsw.edu.au>.
#' @section References:
#' Dunn, P.K., & Smyth, G.K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics 5, 236-244.
#' 
#' Popovic, G. C., Hui, F. K., & Warton, D. I. (2018). A general algorithm for covariance modeling of discrete data. Journal of Multivariate Analysis, 165, 86-100.
#' @section See also:
#' \code{\link{plot.cord}}
#' @import mvabund
#' @export 
#' @examples
#' data(spider)
#' abund <- mvabund(spider$abund)
#' X <- spider$x
#' spider_mod=manyglm(abund~1)
#' spid_lv=cord(spider_mod)
#' plot(spid_lv,biplot = TRUE)
cord <- function(obj, nlv = 2, n.samp = 500, seed = 1) {
    
    if (!is.numeric(seed)) 
        stop("seed must be numeric")
    
    if (floor(n.samp) != ceiling(n.samp)) 
        stop("n.samp must be an integer")
    
    if (floor(nlv) != ceiling(nlv)) 
        stop("nlv must be an integer")
    
    # always same result unless specified otherwise
    set.seed(seed)
    
    # simulate full set of residulas n.samp times
    res = simulate.res.S(obj, n.res = n.samp)
    
    # carry out factor analysis
    out = fact.many(obj, nlv, res)
    
    class(out) = "cord"
    return(out)
    
}
