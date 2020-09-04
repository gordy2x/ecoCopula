#'Simulates new data from a given cord object
#'
#' @param object is a cord object, e.g. from output of \code{cord}
#' @param nsim Number of simulations, defaults to 1.
#' @param seed Random number seed, defaults to a random seed number.
#' @param newdata A data frame in which to look for X covariates with which to simulate.
#' Defaults to the X covariates in the fitted model.
#' @param reshape2matrix logical. Returns a matrix of simulated data if TRUE or nsim == 1.
#' If TRUE and nsim > 1, the simulated data will be appended and returned as a matrix.
#' Returns a 3-d array of simulated data if FALSE and nsim > 1. Each simulation is returned
#' along the third dimension. Defaults to TRUE.
#' @export
#' @examples
#' data(spider)
#' abund <- mvabund(spider$abund)
#' X <- spider$x
#' spider_mod=manyglm(abund~1)
#' spid_lv=cord(spider_mod)
#' simulate(spid_lv)
#' Xnew = data.frame(X[1:10,])
#' simulate(spid_lv, newdata = Xnew)
#' simulate(spid_lv, nsim=2, newdata = Xnew)


simulate.cord = function(object, nsim=1, seed=NULL, 
  newdata=object$obj$data, reshape2matrix=TRUE) {
  
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

  prs = predict_responses(object, newdata)
  newY = simulate_newY(object, nsim, prs)
  newY = reshape_newY(object, nsim, reshape2matrix, newY)

  return (newY)
}

predict_responses = function(object, newdata) {
  if (object$obj$call[[1]] == "manyglm") {
    design.matrix = model.matrix(object$obj$formula[-2], data=newdata)
    coeffs = coef(object$obj)
    
    if (object$obj$family == "negative.binomial") {
      prs = MASS::negative.binomial(theta=1)$linkinv(design.matrix%*%coeffs)
    } else if (object$obj$family == "poisson") {
      prs = poisson(link="log")$linkinv(design.matrix%*%coeffs)
    } else if (object$obj$family == "binomial(link=logit)") {
      prs = binomial(link="logit")$linkinv(design.matrix%*%coeffs)
    } else if (object$obj$family == "binomial(link=cloglog)") {
      prs = binomial(link="cloglog")$linkinv(design.matrix%*%coeffs)
    } else {
      stop("'family'", object$family, "not recognized")
    }
  } else if (object$obj$call[[1]] == "manylm") {
    prs = predict.manylm(object$obj, type = "response")
  } else {
    stop("'class'", class(object$obj), "not supported")
  }
 
  return (prs)
}

simulate_newY = function(object, nsim, prs) {
  nRow = nrow(prs)
  nVar = ncol(prs)
  newY = array(NA, c(nRow, nVar, nsim))

  # simulate multivariate normal random variables
  sim = MASS::mvrnorm(nRow, mu = rep(0, times = nVar), object$sigma)

  # turn simulated variables into abundances
  for (k in 1:nsim) {
    if (object$obj$call[[1]] == "manyglm") {
      for (j in 1:nVar) {
        if (object$obj$family == "negative.binomial") {
          size = object$obj$theta
          newY[,j,k] = qnbinom(pnorm(sim[,j]), size = size[j], mu = prs[,j])
        } else if (object$obj$family == "poisson") {
          newY[,j,k] = qpois(pnorm(sim[,j]), lambda = prs[,j])
        } else if (object$obj$call$family == "binomial") {
          newY[,j,k] = qbinom(pnorm(sim[,j]), size = 1, prob = prs[,j])
        } else {
          stop("'family'", object$family, "not recognized")
        }
      }
    }

    if (object$obj$call[[1]] == "manylm") {
      df.residual = object$obj$df.residual
      sigma2 = apply((object$obj$y - object$obj$fitted)^2, 2, sum)/df.residual
      for (j in 1:nVar) {
        newY[,j,k] = sim[,j] * sqrt(sigma2[j])
      }
      newY[,,k] = newY[,,k] + prs
    }
  }

  return (newY)
}

reshape_newY = function(object, nsim, reshape2matrix, newY) {
  if (reshape2matrix == TRUE || nsim == 1) {
    ylist = lapply(seq(dim(newY)[3]), function(k) newY[,,k])
    ymat = matrix(NA, dim(newY)[1]*nsim, dim(newY)[2])
    ymat = do.call("rbind", ylist)
    out = ymat
  }

  if (reshape2matrix == FALSE && nsim > 1) {
    out = newY
  }

  colnames(out) = colnames(object$obj$y)
  return (out)
}