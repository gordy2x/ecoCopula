
# to import a function use this but with the roxygen comment type #' @importFrom mvabund manyglm

# mid level functions

glasso_opt = function(rho, S.list, full = FALSE, quick = FALSE, start = "cold", cov.last = NULL, prec.last = NULL) {
    P = dim(S.list[[1]])[1]
    J = length(S.list)
    eps = 1e-3
    maxit = 10
    array.S = array(unlist(S.list), c(P, P, J))
    S.init = cov2cor(apply(array.S, c(1, 2), mean))
    weights = rep(1, J)/J
    
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
    
    if (!(quick)) {
        count = 1
        diff = eps + 1
        diff_old=diff+1
        while ((diff > eps & count < maxit) & !any(is.na(Theta.gl[[count]]))) {
            weights = plyr::laply(S.list, L.icov.prop, Theta = Theta.gl[[count]])
            weights = weights/sum(weights)
            ESS=sum(weights)^2/sum(weights^2)
            
            if(ESS<(J/20)|(diff_old<diff))
              stop("this is bad, probable divergence, try removing some columns, increasing n.samp and changing the seed")            # plot(weights)
            diff_old=diff
            count = count + 1
            Sigma.gl[[count]] = cov2cor(apply(array.S, c(1, 2), weighted.mean, w = weights))
            
                A = suppressWarnings(glasso::glasso(Sigma.gl[[count]], rho = rho, penalize.diagonal = FALSE, thr = 1e-08, 
                  start = "warm", w.init = cov.last, wi.init = prec.last))
                cov.last = A$w
                prec.last = A$wi

            
            
            Sigma.gl[[count]] = A$w
            Theta.gl[[count]] = A$wi
            if (!any(is.na(Theta.gl[[count]]))) {
                diff = mean(((Theta.gl[[count]] - Theta.gl[[count - 1]])^2))
            } else {
                diff = mean(((Sigma.gl[[count]] - Sigma.gl[[count - 1]])^2))
            }
            # print(diff)
        }
        
    }
    # print(rho)
    if (full) {
        return(A)
    } else {
        return(A$wi)
    }
}







full.graph.many <- function(manyglm.obj, lambdas, res) {
    S.list = res$S.list
    res = res$res
    P = dim(S.list[[1]])[1]
    N = dim(manyglm.obj$fitted.values)[1]
    
    Th.out = Sig.out = list()
    n.lam = length(lambdas)
    A = glasso_opt(lambdas[1], S.list, full = TRUE, quick = FALSE)
    cov.last = A$w
    prec.last = A$wi
    Th.out[[1]] = prec.last
    Sig.out[[1]] = cov.last
    
    if (n.lam > 1) {
        
        for (i.lam in 2:n.lam) {
            A = glasso_opt(lambdas[i.lam], S.list, full = TRUE, quick = FALSE, start = "warm", cov.last = cov.last, 
                prec.last = prec.last)
            cov.last = A$w
            prec.last = A$wi
            Th.out[[i.lam]] = prec.last
            Sig.out[[i.lam]] = cov.last
            
        }
    }
    
    k = rep(NA, length(Th.out))
    for (j in 1:length(Th.out)) {
        k[j] = (sum(Th.out[[j]] != 0) - P)/2
    }
    k.frac = k/(P * (P - 1)/2)
    BIC.out = logL = NULL
    
    
    logL = plyr::laply(Th.out, ll.icov.all, S.list = S.list, n = N)
    BIC.out = k * log(N) - 2 * logL  #-sum(manyglm.obj$two.loglike)
    AIC.out = k * 2 - 2 * logL  #-sum(manyglm.obj$two.loglike)
    
    
    return(list(k.frac = k.frac, BIC = BIC.out, AIC = AIC.out,logL = logL, Th.out = Th.out, Sig.out = Sig.out))
    
}


fact.many <- function(manyglm.obj, nlv, res) {
    S.list = res$S.list
    res = res$res
    P = dim(S.list[[1]])[1]
    N = dim(manyglm.obj$fitted.values)[1]
    
    
    n.nlv = length(nlv)
    Th.out = Sig.out = Loadings = Scores = list(n.nlv)
    
    for (i.nlv in 1:n.nlv) {
        A = factor_opt(nlv[i.nlv], S.list, full = TRUE, quick = FALSE, N = N)
        Th.out[[i.nlv]] = A$theta
        Sig.out[[i.nlv]] = A$sigma
        Loadings[[i.nlv]] = A$loadings
        res.mean <-  plyr::aaply(plyr::laply(res,function(x) x),c(2,3),weighted.mean,weighs=A$weights)
        # res.mode <- plyr::laply(res,function(x) x)[which(A$weights==max(A$weights)),,]
        Scores[[i.nlv]] = t(as.matrix(A$loadings)) %*% A$theta %*% t(res.mean)
    }
    
    BIC.out = logL = NULL
    k = P * nlv + P - nlv * (nlv - 1)/2
    logL = plyr::laply(Th.out, ll.icov.all, S.list = S.list, n = N)
    BIC.out = k * log(N) - 2 * logL  -sum(manyglm.obj$two.loglike)
    
    
    return(list( loadings = Loadings, scores = Scores, 
                 theta = Th.out, sigma = Sig.out,
                 BIC = BIC.out, logL = logL,
                 obj=manyglm.obj))
    
}

# nlv=2 S.list=list(cor(spider$abund)) N=100
factor_opt = function(nlv, S.list, full = FALSE, quick = FALSE, N) {
    P = dim(S.list[[1]])[1]
    J = length(S.list)
    eps = 1e-10
    maxit = 10
    array.S = array(unlist(S.list), c(P, P, J))
    S.init = cov2cor(apply(array.S, c(1, 2), mean))
    weights = rep(1, J)/J
    
    
    A = factanal(NA, nlv, covmat = S.init, nobs = N, nstart = 3)
    L = A$loadings
    Fac = A$factors
    
    #### need to calculate scores
    
    
    
    Psi = diag(diag(S.init - L %*% t(L)))
    PsiInv = diag(1/diag(Psi))
    Sest = L %*% t(L) + Psi
    Test = solve(Sest)
    
    
    Sigma.gl = Theta.gl = list()
    Sigma.gl[[1]] = Sest
    Theta.gl[[1]] = Test
    A$sigma = Sest
    A$theta = Test
    
    
    if (!(quick)) {
        count = 1
        diff = eps + 1
        while ((diff > eps & count < maxit) & any(!is.na(Theta.gl[[count]]))) {
            weights = plyr::laply(S.list, L.icov.prop, Theta = Theta.gl[[count]])
            weights = weights/sum(weights)
            count = count + 1
            Sigma.gl[[count]] = cov2cor(apply(array.S, c(1, 2), weighted.mean, w = weights))
            
            
            A = factanal(NA, nlv, covmat = Sigma.gl[[count]], nobs = N, nstart = 3)
            L = A$loadings
            Fac = A$factors
            
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



