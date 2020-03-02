simulate.cord = function(fit, nsim = 1, seed = NULL) {
  
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
  
  # generate new data
  ys <- array(NaN, c(nrow(fit$obj$y), ncol(fit$obj$y), nsim))
  if (fit$obj$call[[1]] == "manyglm") {
    prs <- suppressWarnings(predict.manyglm(fit$obj, type = "response")) # warning for family=poisson suppressed
  }
  if (fit$obj$call[[1]] == "manylm") {
    prs <- predict.manylm(fit$obj, type = "response")
  }

  for (k in c(1:nsim)) {
    # simulate multivariate normal random variables
    sim <- mvrnorm(nrow(fit$obj$y), mu = rep(0, times = ncol(fit$obj$y)), fit$sigma[[1]])
    # turn simulated variables into abundances
    if (fit$obj$call[[1]] == "manyglm") {
      if (fit$obj$family == "negative.binomial") {
        size <- fit$obj$theta
        for (j in c(1:ncol(fit$obj$y))) {
          ys[,j,k] <- qnbinom(pnorm(sim[,j]), size = size[j], mu = prs[,j])
        }
      } else if (fit$obj$family == "poisson") {
        for (j in c(1:ncol(fit$obj$y))) {
          ys[,j,k] <- qpois(pnorm(sim[,j]), lambda = prs[,j])
        }
      } else if (fit$obj$call$family == "binomial") {
        for (j in c(1:ncol(fit$obj$y))) {
          ys[,j,k] <- qbinom(pnorm(sim[,j]), size = 1, prob = prs[,j])
        }
      } else {
        stop (paste("'family'", fit$obj$family, "is not supported currently."))
      }
    } else if (fit$obj$call[[1]] == "manylm") {
      df.residual <- fit$obj$df.residual
      sigma2 <- apply((fit$obj$y - fit$obj$fitted)^2, 2, sum)/df.residual
      for (j in c(1:ncol(fit$obj$y))) {
        ys[,j,k] <- qnorm(pnorm(sim[,j]), mean = prs[,j], sd = sqrt(sigma2[j]))
        # ys[,j,k] <- rnorm(ys[,j,k], mean = prs[,j], sd = sqrt(sigma2[j]))
      }
    } else {
      stop (paste("'function'", fit$obj$call[[1]], "is not supported currently."))
    }
  }

  colnames(ys) <- colnames(fit$obj$y)
  return (ys)
}