
#' Fitting Gaussian copula graphical lasso to co-occurence data
#'
#' \code{saint} is used to fit a Gaussian copula graphical model to 
#' multivatiate discrete data, like species co-occurence data in ecology. 
#' This function fits the model and estimates the shrinkage paramater
#' using BIC. Use \code{\link{plot.saint}} to plot the resulting graph.
#'
#' @param obj object of either class \code{\link[mvabund]{manyglm}}, 
#' or  \code{\link[mvabund]{manyany}} with ordinal models \code{\link[ordinal]{clm}}
#' @param lambda vector, values of shrinkage parameter lambda for model 
#' selection (optional, see detail)
#' @param n.lambda integer, number of lambda values 
#' for model selection (optional), ignored if lambda suplied
#' @param n.samp integer, number of sets residuals used for importance sampling 
#' (optional, see detail)
#' @param seed integer, seed for random number generation
#' @section Details:
#' Blah
#' @return Three objects are returned; \code{best_graph} is a list with 
#' paramaters for the 'best' graphical model, chosen by BIC; 
#' \code{all_graphs} is a list with likelihood and BIC for all models along lambda path; 
#' \code{obj} is the input object.
#' @section Author(s):
#' jh
#' @section References:
#' jhg
#' @section See also:
#' \code{\link{plot.saint}}
#' @examples
#' data(spider)
#' abund <- mvabund(spider$abund)
#' X <- spider$x
#' spider.mod=manyglm(abund~X)
#' spid.graph=saint(spider.mod)
#' plot(spid.graph)
#' 
#' 
#' @export
saint <- function(obj, lambda = NULL, n.lambda = 100, n.samp = 500, seed = 1) {
    
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
    logL = ag$logL
    
    
    best = min(which(BIC.graph == min(BIC.graph)))
    Th.best = ag$Th.out[[best]]
    Sig.best = ag$Sig.out[[best]]
    
    graph.out = as.matrix((Th.best != 0) * 1)
    best.graph = list(graph = graph.out, prec = Th.best, cov = Sig.best, Y = obj$y, logL = logL[[best]], 
        sparsity = k.frac[best])
    all.graphs = list(lambda.opt = lambda[best], logL = logL, BIC = BIC.graph, lambda = lambda, k.frac = k.frac)
    out = list(best_graph = best.graph, all_graphs = all.graphs, obj = obj)
    class(out) = "saint"
    return(out)
    
}
