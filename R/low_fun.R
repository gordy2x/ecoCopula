




# low level functions

simulate.res.S <- function(obj, n.res = 200) {
    # pars = params(obj) many.pars = list(pars)[rep(1, n.res)] res = plyr::llply(many.pars, residuals)
    
    many.obj = list(obj)[rep(1, n.res)]
    res = plyr::llply(many.obj, residuals)
    
    if (min(res[[1]]) > -1e-05) {
        # {residuals function for manyany currently output on uniform scale
        res = plyr::llply(res, qnorm)
    }
    S.list = plyr::llply(res, function(x) cov2cor(cov0(x)))  #cov0 assumes 0 mean
    return(list(res = res, S.list = S.list))
}


cov0 = function(Y) {
    n = dim(Y)[1]
    (t(Y) %*% Y)/(n - 1)
}


ll.icov.all <- function(theta, S.list, n) {
    P = dim(S.list[[1]])[1]
    K = length(S.list)
    inside = 0
    for (i in 1:K) {
        scatter = (n - 1) * S.list[[i]]
        inside = inside + exp(-sum(diag(scatter %*% (theta - diag(P))))/2)
    }
    
    ll = -n/2 * log(2 * pi) + n/2 * log(det(theta)) + log(inside/n)
    # notice theta is sigma inverse hence the first minus is a plus
}


print.saint = function(obj) {
    print(paste0("Species interactions from copula model", obj$formula))
    print("Precision matrix")
    print(obj$best.graph$prec)
}

L.icov.prop = function(S, Theta) {
    # The trace of a product can be rewritten as the sum of entry-wise products of elements
    exp(-1/2 * sum(S * (Theta - diag(dim(Theta)[1]))))
}


params.manyglm <- function(object) {
    
    # code pilfered from residuals.manyglm
    n.rows = NROW(object$y)
    n.vars = NCOL(object$y)
    params = list()
    if (object$family == "negative.binomial") {
        pfn = "pnbinom"
        for (i.var in 1:n.vars) params[[i.var]] = list(q = object$y[, i.var], mu = object$fitted[, 
            i.var], size = object$theta[i.var])
    } else if (object$family == "poisson") {
        pfn = "ppois"
        for (i.var in 1:n.vars) params[[i.var]] = list(q = object$y[, i.var], lambda = object$fitted[, 
            i.var])
    } else if (substr(object$family, 1, 3) == "bin" || "clo") {
        pfn = "pbinom"
        for (i.var in 1:n.vars) params[[i.var]] = list(q = object$y[, i.var], size = 1, prob = object$fitted[, 
            i.var])
    } else if (object$family == "gaussian") {
        pfn = "pnorm"
        df.residual = n.rows - dim(coef(object))[1]
        sigma2 = apply((object$y - object$fitted)^2, 2, sum)/df.residual
        for (i.var in 1:n.vars) params[[i.var]] = list(q = object$y[, i.var], mu = object$fitted[, 
            i.var], sd = sqrt(sigma2[i.var]))
    } else stop(paste("'family'", object$family, "not recognized"))
    
    qupper = qlower = matrix(NA, n.rows, n.vars)
    for (i.var in 1:n.vars) {
        param.minus = params[[i.var]]
        param.minus$q = params[[i.var]]$q - 1e-06
        qupper[, i.var] = do.call(pfn, params[[i.var]])
        qlower[, i.var] = do.call(pfn, param.minus)
    }
    out = list(qlower = qlower, qupper = qupper)
    class(out) = "params"
    return(out)
}


params.manyany <- function(object) {
    # code pilfered from residuals.manyany
    
    params = object$params
    n.rows = length(params[[1]]$q)
    n.vars = length(params)
    if (length(object$family) == 1) 
        family = rep(object$family, n.vars) else family = object$family
    qupper = qlower = matrix(NA, n.rows, n.vars)
    for (i.var in 1:n.vars) {
        if (family[[i.var]]$family == "ordinal") {
            qupper[, i.var] = params[[i.var]]$mu$cprob2
            qlower[, i.var] = params[[i.var]]$mu$cprob1
        } else stop(paste("'family'", object$family, "not recognized"))
    }
    out = list(qlower = qlower, qupper = qupper)
    class(out) = "params"
    return(out)
}

params <- function(x, ...) {
    UseMethod("params", x)
}

residuals.params <- function(params) {
    # code pilfered from residuals.manyglm
    tol = 1e-08
    qupper = params$qupper
    qlower = params$qlower
    n.rows = NROW(qlower)
    n.vars = NCOL(qlower)
    u = matrix(runif(n.rows * n.vars), n.rows, n.vars)
    resids = u * pmax(tol^3, qupper) + (1 - u) * pmin(1 - tol^3, qlower)
    return(qnorm(resids))
}

