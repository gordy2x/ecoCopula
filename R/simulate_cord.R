#'Simulates new data from a given cord object
#'
#' @param object is a cord object, e.g. from output of \code{cord}
#' @param nsim Number of simulations, defaults to 1. If nsim > 1, the simulated data will be
#' appended.
#' @param seed Random number seed, defaults to a random seed number.
#' @param newdata A data frame in which to look for X covariates with which to simulate.
#' Defaults to the X covariates in the fitted model.
#' @export
#' @examples
#' data(spider)
#' abund = mvabund(spider$abund)
#' X = data.frame(spider$x)
#'
#' spider_mod = manyglm(abund~1)
#' spid_lv = cord(spider_mod)
#' simulate(spid_lv)
#'
#' spider_mod_X = manyglm(abund ~ soil.dry + bare.sand, data=X)
#' spid_lv_X = cord(spider_mod_X)
#' Xnew = X[1:10,]
#' simulate(spid_lv_X, newdata = Xnew)
#' simulate(spid_lv_X, nsim=2, newdata = Xnew)


simulate.cord = function(object, nsim=1, seed=NULL, newdata=object$obj$data) {
  
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

  newdata = reshape_newdata(object, nsim, newdata)
  prs = suppressWarnings(
      predict(object$obj, type = "response", newdata = newdata)
  ) # warning for family=poisson suppressed
  newY = simulate_newY(object, prs)

  return (newY)
}

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
        newdata = data.frame(rep(newdata, nsim))

        if (ncol(newdata) > 1) {
          newdata = reshape(newdata, direction = "long", varying=1:ncol(newdata))
        }
      }
    }
  }

 return (newdata)
}

simulate_newY = function(object, prs) {
  nRow = nrow(prs)
  nVar = ncol(prs)
  newY = matrix(NA, nRow, nVar)

  # simulate multivariate normal random variables
  sim = MASS::mvrnorm(nRow, mu = rep(0, times = nVar), object$sigma)

  # turn simulated variables into abundances
  if (object$obj$call[[1]] == "manyglm") {
    for (iVar in 1:nVar) {
      if (object$obj$family == "negative.binomial") {
        size = object$obj$theta
        newY[,iVar] = qnbinom(pnorm(sim[,iVar]), size = size[iVar], mu = prs[,iVar])
      } else if (object$obj$family == "poisson") {
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
