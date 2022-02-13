#' Fitting Gaussian copula graphical lasso to co-occurrence data
#'
#' \code{cgr} is used to fit a Gaussian copula graphical model to 
#' multivariate discrete data, like species co-occurrence data in ecology. 
#' This function fits the model and estimates the shrinkage parameter
#' using BIC. Use \code{\link{plot.cgr}} to plot the resulting graph.
#'
#' @param obj object of either class \code{\link[mvabund]{manyglm}}, 
#' or  \code{\link[mvabund]{manyany}} with ordinal models \code{\link[ordinal]{clm}}
#' @param lambda vector, values of shrinkage parameter lambda for model 
#' selection (optional, see detail)
#' @param n.lambda integer, number of lambda values 
#' for model selection (default = 100), ignored if lambda supplied
#' @param n.samp integer (default = 500), number of sets residuals used for importance sampling 
#' (optional, see detail)
#' @param method method for selecting shrinkage parameter lambda, either "BIC" (default) or "AIC"
#' @param seed integer (default = 1), seed for random number generation (optional, see detail)
#' @section Details:
#' \code{cgr} is used to fit a Gaussian copula graphical model to multivariate discrete data, such as co-occurrence (multi species) data in ecology. The model is estimated using importance sampling with \code{n.samp} sets of randomised quantile or "Dunn-Smyth" residuals (Dunn & Smyth 1996), and the \code{\link{glasso}} package for fitting Gaussian graphical models. Models are fit for a path of values of the shrinkage parameter \code{lambda} chosen so that both completely dense and sparse models are fit. The \code{lambda} value for the \code{best_graph} is chosen by BIC (default) or AIC.  The seed is controlled so that models with the same data and different predictors can be compared.  
#' @return Three objects are returned; 
#' \code{best_graph} is a list with parameters for the 'best' graphical model, chosen by the chosen \code{method}; 
#' \code{all_graphs} is a list with likelihood, BIC and AIC for all models along lambda path; 
#' \code{obj} is the input object.
#' @section Author(s):
#' Gordana Popovic <g.popovic@unsw.edu.au>.
#' @section References:
#' Dunn, P.K., & Smyth, G.K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics 5, 236-244.
#' 
#' Popovic, G. C., Hui, F. K., & Warton, D. I. (2018). A general algorithm for covariance modeling of discrete data. Journal of Multivariate Analysis, 165, 86-100.
#' @section See also:
#' \code{\link{plot.cgr}}
#' @examples
#' abund <- spider$abund[,1:5]
#' spider_mod <- stackedsdm(abund,~1, data = spider$x, ncores=2) 
#' spid_graph=cgr(spider_mod)
#' plot(spid_graph,pad=1)
#' @import mvabund
#' @export 
cgr <- function(obj, lambda = NULL, n.lambda = 100, 
                  n.samp = 500, method="BIC", seed = NULL) {
    
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
    
    if (!class(obj)[1] %in% c("manyany", "manyglm","manylm","stackedsdm")) 
        stop("please supply an manyglm, manylm, manyany or stackedsdm object")
    
    if (class(obj)[1] == "manyany" & class(obj)[2] != "clm")
        warning("cgr function is only tested on manyany with clm or tweedie")
    
    # always same result unless specified otherwise
    set.seed(seed)
    
    # simulate full set of residuals n.samp times
    res = simulate.res.S(obj, n.res = n.samp)
    
    #this chunk of code finds a path of lambda values to 
    #explore graphs ranging from completely sparse to completely dense.
    if (is.null(lambda)) {
        # starting values for log10 lambda
        current = seq(-6, 2, length.out = 10)
        
        # proportion of non-zero cond indep
        sparse_fits = full.graph.many(obj, c(0,10^current), res)
        k.frac=sparse_fits$k.frac
        
        # which lambdas give full and empty matrices
        new.min = max(which(k.frac == 1))
        new.max = min(which(k.frac == 0))
        
        # now a small sequence of lambda values, between the two
        lambda = c(0,10^seq(current[new.min], current[new.max], length.out = 10))
        sparse_fits = full.graph.many(obj, lambda, res)
        
        
        #then find the sweet spot by removing the 5 largest of the 10 values above
        if(method=="BIC"){
          sub_lam=lambda[sparse_fits$BIC<kth_largest(sparse_fits$BIC,5)]
        }else if (method=="AIC"){
          sub_lam=lambda[sparse_fits$AIC<kth_largest(sparse_fits$AIC,5)]
        }else{
          stop("lambda selection method can only be \"AIC\" or \"BIC\" ")
        }
        
        # now a larger sequence but just between those
        min_lam=min(sub_lam)
        if(min_lam>0){
          lambda = 10^c( seq(log10(min_lam), log10(max(sub_lam)), length.out = n.lambda))
        } else{
          min_lam=min(sub_lam[sub_lam>0])
          lambda = c(0,10^c( seq(log10(min_lam), log10(max(sub_lam)), length.out = (n.lambda-1))))
        }

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
    
    
    # find raw correlation matrix by averaging over unweighted residuals
    P = dim(res$S.list[[1]])[1]
    array.S = array(unlist(res$S.list), c(P, P, n.samp))
    precov = cov2cor(apply(array.S, c(1, 2), mean))
    colnames(precov)=rownames(precov)=colnames(obj$y)
    
    if(any(class(obj) == "manyany")){
      labs <- names(obj$params)
    }else{
      labs <- colnames(obj$y)
    }
    
    Th.best = ag$Th.out[[best]]
    Sig.best = ag$Sig.out[[best]]
    colnames(Th.best)=rownames(Th.best)=colnames(Sig.best)=rownames(Sig.best)=labs
    part_cor = -cov2cor(Th.best)
    g<-graph_from_partial(part_cor)
    #outputs
    
    graph.out = as.matrix((Th.best != 0) * 1)
    best.graph = list(graph = graph.out, prec = Th.best, cov = Sig.best, part=part_cor, Y = obj$y, logL = logL[[best]], 
        sparsity = k.frac[best],igraph_out=g)
    all.graphs = list(lambda.opt = lambda[best], logL = logL, BIC = BIC.graph, 
                      AIC = AIC.graph, lambda = lambda, k.frac = k.frac)
    raw <- list(cov=precov)
    out = list(best_graph = best.graph,raw=raw, all_graphs = all.graphs, obj = obj)
    class(out) = "cgr"
    return(out)
    
}

#' glasso path
#'
#' Fits IRW glasso path
#'
#' @param obj fitted manyglm or stackedsdm object
#' @param lambdas vector of shrinkage parameters
#' @param res object containing residuals and their empirical covariance matrices
#'
#' @noRd
full.graph.many <- function(obj, lambdas, res) {
  
  S.list = res$S.list
  res = res$res
  P = dim(S.list[[1]])[1]
  N = dim(obj$fitted)[1]
  
  Th.out = Sig.out = list()
  n.lam = length(lambdas)
  
  #first lambda fit, cold start
  A = glasso_opt(lambdas[1], S.list, full = TRUE, quick = FALSE)
  cov.last = A$w
  prec.last = A$wi
  Th.out[[1]] = prec.last
  Sig.out[[1]] = cov.last
  
  if (n.lam > 1) {
    
    for (i.lam in 2:n.lam) {
      
      #warm start for other lambdas
      A = glasso_opt(lambdas[i.lam], S.list, full = TRUE, quick = FALSE, start = "warm", cov.last = cov.last, 
                     prec.last = prec.last)
      cov.last = A$w
      prec.last = A$wi
      Th.out[[i.lam]] = prec.last
      Sig.out[[i.lam]] = cov.last
      
    }
  }
  
  # calculate likelihood and AIC/BIC
  k = rep(NA, length(Th.out))
  for (j in 1:length(Th.out)) {
    k[j] = (sum(Th.out[[j]] != 0) - P)/2
  }
  k.frac = k/(P * (P - 1)/2)
  BIC.out = logL = NULL
  
  
  logL = plyr::laply(Th.out, ll.icov.all, S.list = S.list, n = N)
  BIC.out = k * log(N) - 2 * logL  #-sum(obj$two.loglike)
  AIC.out = k * 2 - 2 * logL  #-sum(obj$two.loglike)
  
  
  return(list(k.frac = k.frac, BIC = BIC.out, AIC = AIC.out,logL = logL, Th.out = Th.out, Sig.out = Sig.out))
  
}

#' glasso iteration
#'
#' Fits graphical lasso for single shrinkage parameter value
#'
#' @param rho shrinkage parameter
#' @param S.list list of empirical covariance matrices for sets of residuals
#' @param full boolean, should whole fit be returned or just precision matrix
#' @param quick boolean, if true only one iteration is done
#' @param start should glasso use warm start. either "cold", or "warm"
#' @param cov.last warm start covariance estimate
#' @param prec.last warm start precision estimate
#'
#' @noRd
glasso_opt = function(rho, S.list, full = FALSE, quick = FALSE, start = "cold", cov.last = NULL, prec.last = NULL) {
  
  P = dim(S.list[[1]])[1]
  J = length(S.list)
  eps = 1e-3
  maxit = 10
  
  #list of correlation of each set of residuals
  array.S = array(unlist(S.list), c(P, P, J)) 
  
  #initialize covariance and weights
  S.init = cov2cor(apply(array.S, c(1, 2), mean)) 
  weights = rep(1, J)/J 
  
  #starting parameters if warm start
  if (start == "warm") {
    A = suppressWarnings(glasso::glasso(S.init, rho = rho, penalize.diagonal = FALSE, thr = 1e-08, start = "warm", 
                                        w.init = cov.last, wi.init = prec.last))
  } else {
    A = suppressWarnings(glasso::glasso(S.init, rho = rho, penalize.diagonal = FALSE, thr = 1e-08))
  }
  
  cov.last = A$w
  prec.last = A$wi
  
  Sigma.gl = Theta.gl = list()
  Sigma.gl[[1]] = A$w
  Theta.gl[[1]] = A$wi
  
  #if quick only one iteration, otherwise iterate
  if (!(quick)) {
    count = 1
    diff = eps + 1
    diff_old=diff+1
    while ((diff > eps & count < maxit) & !any(is.na(Theta.gl[[count]]))) {
      
      #recalculate weights
      weights = plyr::laply(S.list, L.icov.prop, Theta = Theta.gl[[count]])
      weights = weights/sum(weights)
      
      #effective sample size check
      ESS=sum(weights)^2/sum(weights^2)
      if(ESS<(J/20)|(diff_old<diff))
        stop("this is bad, probable divergence, try removing some columns, increasing n.samp and changing the seed")
      
      diff_old=diff
      count = count + 1
      
      #recalculate weighted covariance matrix
      Sigma.gl[[count]] = cov2cor(apply(array.S, c(1, 2), weighted.mean, w = weights))
      
      #fit glasso and extract estimates
      A = suppressWarnings(glasso::glasso(Sigma.gl[[count]], rho = rho, penalize.diagonal = FALSE, thr = 1e-08, 
                                          start = "warm", w.init = cov.last, wi.init = prec.last))
      cov.last = A$w
      prec.last = A$wi
      
      #convergence
      Sigma.gl[[count]] = A$w
      Theta.gl[[count]] = A$wi
      if (!any(is.na(Theta.gl[[count]]))) {
        diff = mean(((Theta.gl[[count]] - Theta.gl[[count - 1]])^2))
      } else {
        diff = mean(((Sigma.gl[[count]] - Sigma.gl[[count - 1]])^2))
      }
    }
    
  }
  if (full) {
    return(A)
  } else {
    return(A$wi)
  }
}




