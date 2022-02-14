
utils::globalVariables(c("cov2cor", "factanal", "weighted.mean","coef","plot",
                         "abline","text","arrows","runif","qnorm","residuals"))


#' Simulate residuals
#'
#' Simulates multiple residuals
#'
#' @param obj manyglm or stackesdm object
#' @param n.res number of residual sets to simulate
#'
#' @noRd
simulate.res.S <- function(obj, n.res = 200) {

    many.obj = list(obj)[rep(1, n.res)]
    res = plyr::llply(many.obj, residuals)
    res = plyr::llply(res, fix_inf)
    
    # residuals function for manyany currently output on uniform scale
    if (min(res[[1]]) > -1e-05) {
        res = plyr::llply(res, qnorm)
    }
    
    S.list = plyr::llply(res, function(x) cov2cor(cov0(x)))  #cov0 assumes 0 mean
    return(list(res = res, S.list = S.list))
}

#' Cov assuming 0 mean
#'
#' empirical covariance assuming mean 0
#'
#' @param Y matrix
#'
#' @noRd
cov0 = function(Y) {
    n = dim(Y)[1]
    (t(Y) %*% Y)/(n - 1)
}

#' multiple log likelihood 
#'
#' Calculate log likelihood with multiple empirical covariance matrices and one estimated precision
#'
#' @param theta estimated precision
#' @param S.list list of empirical covariance matrices for sets of residuals
#' @param nobs sample size
#'
#' @noRd
ll.icov.all <- function(theta, S.list, nobs) {
    P = dim(S.list[[1]])[1]
    K = length(S.list)
    inside = 0
    for (i in 1:K) {
        scatter = (nobs - 1) * S.list[[i]] #scatter matrix
        inside = inside + exp(-sum(diag(scatter %*% (theta - diag(P))))/2)
    }
    
    ll = -nobs/2 * log(2 * pi) + nobs/2 * log(det(theta)) + log(inside/nobs)
    # notice theta is sigma inverse hence the first minus is a plus
}


#' single log likelihood unnormalised
#'
#' Calculate log likelihood (unnormalised) with empirical covariance matrix and estimated precision
#'
#' @param S empirical covariance matrix
#' @param Theta estimated precision
#'
#' @noRd
L.icov.prop = function(S, Theta) {
    
    # The trace of a product can be rewritten as the sum of entry-wise products of elements
    exp(-1/2 * sum(S * (Theta - diag(dim(Theta)[1]))))
}

#' kth larges element
#'
#' Find k'th largest element of a vector
#'
#' @param x vector
#' @param k integer
#'
#' @noRd
kth_largest<-function(x,k=1){
    n <- length(x)
    sort(x,partial=n-k+1)[n-k+1]
}


#' reshape newdata
#'
#' Reshapes newdata
#'
#' @param object 
#' @param nsim 
#' @param newdata
#'
#' @noRd
reshape_newdata = function(object, nsim, newdata) {
    if (formula(object$obj)[-2] == ~1) {
        newdata = data.frame(array(1, nrow(newdata)*nsim))
    } else {
        if (nsim == 1) {
            newdata = newdata
        } else {
            if (ncol(newdata) > 1) {
                newdata = do.call("rbind", replicate(nsim, newdata, simplify = FALSE))
            } else {
                name = colnames(newdata)
                colnames(newdata) = "V1"
                newdata = data.frame(rep(newdata, nsim))
                
                if (ncol(newdata) > 1) {
                    newdata = reshape(newdata, direction = "long", varying=1:ncol(newdata))["V1"]
                }
                colnames(newdata) = name
            }
        }
    }
    
    return (newdata)
}


#' get size
#'
#' Get size 
#'
#' @param object 
#'
#' @noRd
get_size = function(object) {
    if (object$obj$call[[1]] == "manyglm") {
        size = object$obj$theta
    }
    
    if (as.character(object$obj$call[[1]]) == "stackedsdm") {
        fits = object$obj$fits
        size = sapply(fits, function(x) x[['fit']][['theta']])
    }
    
    return (size)
}


#' check family
#'
#' check object family is all the same
#'
#' @param object 
#'
#' @noRd
check_object_family = function(object) {
    if (as.character(object$obj$call[[1]]) == "stackedsdm") {
        if (length(unique(object$obj$family)) > 1) {
            stop("Multiple distributions not currently supported")
        }
    }
}


#' Graph from partial correlation
#'
#' Create igraph object from partial correlation matrix
#'
#' @param partial 
#'
#' @noRd
graph_from_partial<-function(partial){
    diag(partial)<-0
    vertex <- data.frame(name = colnames(partial))
    Edges=igraph::as_data_frame(igraph::graph_from_adjacency_matrix(partial,mode="undirected",weighted = TRUE))
    if(ncol(Edges)==3){
        colnames(Edges)[3]<- "partcor"
    }
    
    igraph::graph_from_data_frame(Edges, directed=FALSE, vertices=vertex)
}

#' Fix infinite residuals
#'
#' Set large residuals to max/min value
#'
#' @param mat a matrix 
#' @param lim the max/min value
#'
#' @noRd
fix_inf<-function(mat, lim=5){
    mat[mat >  lim] = lim
    mat[mat < (-lim)] = -lim
    mat
}

#' get predictions
#'
#' Get predicted responses
#'
#' @param object
#' @param newdata
#'
#' @noRd
predict_responses = function(object, newdata) {
    prs = try(
        suppressWarnings(
            predict(object$obj, type = "response", newdata = newdata)
        ), # warning for family=poisson suppressed
    silent = TRUE)

    if (inherits(prs, "try-error")) {
        coeffs = object$obj$coefficients
        offset = object$obj$offset
        design.matrix = model.matrix(object$obj$formula[-2], data = newdata)
        xbeta = design.matrix %*% coeffs + offset

        if (object$obj$family == "negative.binomial") {
          prs = MASS::negative.binomial(theta = 1)$linkinv(xbeta)
        } else if (object$obj$family == "poisson") {
          prs = poisson(link = "log")$linkinv(xbeta)
        } else if (object$obj$family == "binomial(link=logit)") {
          prs = binomial(link = "logit")$linkinv(xbeta)
        } else if (object$obj$family == "binomial(link=cloglog)") {
          prs = binomial(link = "cloglog")$linkinv(xbeta)
        } else {
          stop("'family'", object$obj$family, "not recognized")
        }
    }

  return (prs)
}
