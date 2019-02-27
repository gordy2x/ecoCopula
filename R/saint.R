#' Fitting Gaussian copula graphical lasso to co-occurence data
#'
#' \code{saint} is used to fit a Gaussian copula graphical model to 
#' multivatiate discrete data, like species co-occurence data in ecology. 
#' This function fits the model and estimates the shrinkage parameter
#' using BIC. Use \code{\link{plot.saint}} to plot the resulting graph.
#'
#' @param obj object of either class \code{\link[mvabund]{manyglm}}, 
#' or  \code{\link[mvabund]{manyany}} with ordinal models \code{\link[ordinal]{clm}}
#' @param lambda vector, values of shrinkage parameter lambda for model 
#' selection (optional, see detail)
#' @param n.lambda integer, number of lambda values 
#' for model selection (default = 100), ignored if lambda supplied
#' @param n.samp integer (default = 500), number of sets residuals used for importance sampling 
#' (optional, see detail)
#' @param seed integer (default = 1), seed for random number generation (optional, see detail)
#' @section Details:
#' \code{saint} is used to fit a Gaussian copula graphical model to multivariate discrete data, such as co-occurence (multi species) data in ecology. The model is estimated using importance sampling with \code{n.samp} sets of randomised quantile or "Dunn-Smyth" residuals (Dunn & Smyth 1996), and the \code{\link{glasso}} package for fitting Gaussian graphical models. Models are fit for a path of values of the shrinkage parameter \code{lambda} chosen so that both completely dense and sparse models are fit. The \code{lambda} value for the \code{best_graph} is chosen by BIC.  The seed is controlled so that models with the same data and different predictors can be compared.  
#' @return Three objects are returned; 
#' \code{best_graph} is a list with parameters for the 'best' graphical model, chosen by BIC; 
#' \code{all_graphs} is a list with likelihood and BIC for all models along lambda path; 
#' \code{obj} is the input object.
#' @section Author(s):
#' Gordana Popovic <g.popovic@unsw.edu.au>.
#' @section References:
#' Dunn, P.K., & Smyth, G.K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics 5, 236-244.
#' 
#' Popovic, G. C., Hui, F. K., & Warton, D. I. (2018). A general algorithm for covariance modeling of discrete data. Journal of Multivariate Analysis, 165, 86-100.
#' @section See also:
#' \code{\link{plot.saint}}
#' @examples
#' data(spider)
#' abund <- mvabund(spider$abund)
#' X <- spider$x
#' spider_mod=manyglm(abund~X)
#' spid_graph=saint(spider_mod)
#' plot(spid_graph,pad=1)
#' 
#' 
#' @export
saint <- function(obj, lambda = NULL, n.lambda = 100, 
                  n.samp = 500, seed = 1, method="BIC") {
    
    if (!is.numeric(seed)) 
        stop("seed must be numeric")
    
    if (floor(n.samp) != ceiling(n.samp)) 
        stop("n.samp must be an integer")
    
    if (floor(n.lambda) != ceiling(n.lambda)) 
        stop("n.lambda must be an integer")
    
    if (!(is.numeric(lambda) | is.null(lambda))) 
        stop("lambda must be numeric")
    
    if (!is.null(lambda)) 
        warning(" 'best' model selected among supplied lambda only")
    
    if (any(lambda < 0)) 
        stop("lambda must be non negative")
    
    if (!class(obj)[1] %in% c("manyany", "manyglm")) 
        stop("please supply an manyglm, or manyany (clm) object")
    
    if (class(obj)[1] == "manyany" & class(obj)[2] != "clm") 
        stop("saint function only supposts manyany with clm")
    
    # always same result unless specified otherwise
    set.seed(seed)
    
    # simulate full set of residulas n.samp times
    res = simulate.res.S(obj, n.res = n.samp)
    
    #this chunk of code finds a path of lambda values to 
    #explore graphs ranging from completely sparse to completely dense.
    if (is.null(lambda)) {
        # starting values for log10 lambda
        current = c(-16, seq(-10, 10, length.out = 20))
        
        # proportion of non-zero cond indep
        k.frac = full.graph.many(obj, 10^current, res)$k.frac
        # which lambdas give full and empty matrices
        new.min = max(which(k.frac == 1))
        new.max = min(which(k.frac == 0))
        
        # now sequence of lambda values, between the two
        lambda = 10^c(-16, seq(current[new.min], current[new.max], length.out = n.lambda))
    } else {
        n.lambda = length(lambda)
    }
    ag = full.graph.many(obj, lambda, res)
    k.frac = ag$k.frac
    BIC.graph = ag$BIC
    AIC.graph = ag$AIC
    logL = ag$logL
    
    #determine best graph by BIC
    if(method=="BIC"){
      best = min(which(BIC.graph == min(BIC.graph)))
    }else if (method=="AIC"){
      best = min(which(AIC.graph == min(AIC.graph)))
    }else{
      stop("lambda selection method can only be \"AIC\" or \"BIC\" ")
    }
    
    
    Th.best = ag$Th.out[[best]]
    Sig.best = ag$Sig.out[[best]]
    
    #outputs
    graph.out = as.matrix((Th.best != 0) * 1)
    best.graph = list(graph = graph.out, prec = Th.best, cov = Sig.best, Y = obj$y, logL = logL[[best]], 
        sparsity = k.frac[best])
    all.graphs = list(lambda.opt = lambda[best], logL = logL, BIC = BIC.graph, AIC = AIC.graph, lambda = lambda, k.frac = k.frac)
    out = list(best_graph = best.graph, all_graphs = all.graphs, obj = obj)
    class(out) = "saint"
    return(out)
    
}
