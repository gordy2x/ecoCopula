#' Model based ordination with Gaussian copulas
#'
#' @param  obj object of either class \code{\link[mvabund]{manyglm}}, 
#' or  \code{\link[mvabund]{manyany}} with ordinal models \code{\link[ordinal]{clm}}
#' @param nlv number of latent variables (default = 2, for plotting on a scatterplot)
#' @param n.samp integer (default = 500), number of sets residuals used for importance sampling 
#' (optional, see detail)
#' @param seed integer (default = 1), seed for random number generation (optional, see detail)

#' @section Details:
#' \code{cord} is used to fit a Gaussian copula factor analytic model to multivariate discrete data, such as co-occurrence (multi species) data in ecology. The model is estimated using importance sampling with \code{n.samp} sets of randomised quantile or "Dunn-Smyth" residuals (Dunn & Smyth 1996), and the \code{\link{factanal}} function. The seed is controlled so that models with the same data and different predictors can be compared.  
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
#' X <- as.data.frame(spider$x)
#' abund <- spider$abund
#' spider_mod <- stackedsdm(abund,~1, data = X) 
#' spid_lv=cord(spider_mod)
#' plot(spid_lv,biplot = TRUE)
cord <- function(obj, nlv = 2, n.samp = 500, seed = NULL) {
    
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
    
    if (length(nlv) != 1) 
        stop("nlv must be an integer")
    
    if (floor(nlv) != ceiling(nlv)) 
        stop("nlv must be an integer")
    
    # always same result unless specified otherwise
    set.seed(seed)
    
    # simulate full set of residuals n.samp times
    res = simulate.res.S(obj, n.res = n.samp)
    
    # carry out factor analysis
    S.list = res$S.list
    res = res$res
    P = dim(S.list[[1]])[1]
    N = dim(obj$fitted)[1]
    A = factor_opt(nlv, S.list, full = TRUE, quick = FALSE, nobs = N)
    
    #extract elements
    Th.out = A$theta
    Sig.out=A$sigma
    colnames(Sig.out)=rownames(Sig.out)=colnames(Th.out)=rownames(Th.out)=colnames(obj$y)
    
    res.mean <-  plyr::aaply(plyr::laply(res,function(x) x),c(2,3),weighted.mean,weighs=A$weights)
    # res.mode <- plyr::laply(res,function(x) x)[which(A$weights==max(A$weights)),,]
    Scores = t(as.matrix(A$loadings)) %*% A$theta %*% t(res.mean)
    
    BIC.out = logL = NULL
    k = P * nlv + P - nlv * (nlv - 1)/2
    logL = ll.icov.all(Th.out, S.list = S.list, nobs = N)
    BIC.out = k * log(N) - 2 * logL  -sum(obj$two.loglike)
    
    out=list( loadings = A$loadings[], scores = t(Scores), 
              sigma = Sig.out, theta=Th.out,
              BIC = BIC.out, logL = logL,
              obj=obj)
    
    
    class(out) = "cord"
    return(out)
    
}

#' IRW factanal
#'
#' Fits iteratively reweighted factanal for one choice of number of latent variables
#'
#' @param nlv number of latent variables
#' @param S.list list of empirical covariance matrices for sets of residuals
#' @param quick boolean, if true only one iteration is done
#' @param nobs number of observations in the data
#'
#' @noRd
factor_opt = function(nlv, S.list, full = FALSE, quick = FALSE, nobs) {
    
    P = dim(S.list[[1]])[1]
    J = length(S.list)
    eps = 1e-10
    maxit = 10
    array.S = array(unlist(S.list), c(P, P, J))
    
    #initialise weighted covariance and weights
    S.init = cov2cor(apply(array.S, c(1, 2), mean))
    weights = rep(1, J)/J
    
    #fit factor analysis
    A = factanal(NA, nlv, covmat = S.init, nobs = nobs, nstart = 3)
    L = A$loadings
    Fac = A$factors
    
    # calculate covariance and precision matrices
    Psi = diag(diag(S.init - L %*% t(L)))
    PsiInv = diag(1/diag(Psi))
    Sest = L %*% t(L) + Psi
    Test = solve(Sest)
    
    
    Sigma.gl = Theta.gl = list()
    Sigma.gl[[1]] = Sest
    Theta.gl[[1]] = Test
    A$sigma = Sest
    A$theta = Test
    
    # if not quick, iterate till convergence
    if (!(quick)) {
        count = 1
        diff = eps + 1
        while ((diff > eps & count < maxit) & any(!is.na(Theta.gl[[count]]))) {
            
            #recalculate weights and weighted covariance
            weights = plyr::laply(S.list, L.icov.prop, Theta = Theta.gl[[count]])
            weights = weights/sum(weights)
            count = count + 1
            Sigma.gl[[count]] = cov2cor(apply(array.S, c(1, 2), weighted.mean, w = weights))
            
            #factor analysis
            A = factanal(NA, nlv, covmat = Sigma.gl[[count]], nobs = nobs, nstart = 3)
            L = A$loadings
            Fac = A$factors
            
            # calculate covariance and precision matrices and save
            Psi = diag(diag(S.init - L %*% t(L)))
            PsiInv = diag(1/diag(Psi))
            Sest = L %*% t(L) + Psi
            Test = solve(Sest)
            
            Sigma.gl = Theta.gl = list()
            Sigma.gl[[count]] = Sest
            Theta.gl[[count]] = Test
            A$sigma = Sest
            A$theta = Test
            A$weights = weights
            
            if (any(!is.na(Theta.gl[[count]]))) {
                diff = sum(((Theta.gl[[count]] - Theta.gl[[count - 1]])^2)/(P^2))
            } else {
                diff = sum(((Sigma.gl[[count]] - Sigma.gl[[count - 1]])^2)/(P^2))
            }
        }
        
    }
    
    
    
    if (full) {
        return(A)
    } else {
        return(A$theta)
    }
}
