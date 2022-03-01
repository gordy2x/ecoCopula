#' Fitted values from a stackedsdm object
#'
#' @param object An object of class \code{stackedsdm}
#' @param ... Not used
#' @section Details:
#' Extracts the fitted values from \code{stackedsdm} object. 
#' @return A matrix of fitted values.
#' @section Author(s):
#' Francis K.C. Hui <francis.hui@anu.edu.au>.
#' @export fitted.stackedsdm
#' @export 
#' @examples
#' library(mvabund)
#' data(spider)
#' X <- spider$x
#' abund <- spider$abund
#'
#' # Example 1: Simple example
#' myfamily <- "negative.binomial"
#' # Example 1: Funkier example where Species are assumed to have different distributions
#' # Fit models including all covariates are linear terms, but exclude for bare sand
#' fit0 <- stackedsdm(abund, formula_X = ~. -bare.sand, data = X, family = myfamily, ncores=2)
#' fitted(fit0)
#'
#' # Example 2: Funkier example where Species are assumed to have different distributions
#' abund[,1:3] <- (abund[,1:3]>0)*1 # First three columns for presence absence
#' myfamily <- c(rep(c("binomial"), 3),
#'               rep(c("negative.binomial"), (ncol(abund)-3)))
#' fit0 <- stackedsdm(abund, formula_X = ~ bare.sand, data = X, family = myfamily, ncores=2)
#' fitted(fit0)
fitted.stackedsdm <- function(object, ...) {
  return(object$fitted)
}


#' Predictions from a stackedsdm object
#'
#' @param object An object of class \code{stackedsdm}
#' @param newdata Optionally, a data frame in which to look for variables with which to predict.  If omitted, the covariates from the existing dataset are used.
#' @param type The type of prediction required.  This can be supplied as either a single character string, when is applied to all species, or a vector of character strings of the same length as \code{ncol(object$y)} specifying the type of predictions desired for each species. The exact type of prediction allowed depends precisely on the distribution, but for many there is at least \code{"link"} which is on the scale of the linear predictors, and \code{"response"} which is on the scale of the response variable. The values of this argument can be abbreviated.
#' @param se.fit Logical switch indicating if standard errors are required.
#' @param na.action Function determining what should be done with missing values in \code{newdata}. The default is to predict \code{NA}..
#' @param ... not used
#' @section Details:
#'  This function simply applies a for loop, cycling through each fitted model from the \code{stackedsdm} object and then attempting to construct the relevant predictions by applying the relevant \code{predict} method. Please keep in mind no formatting is done to the predictions.
#' @return A list where the k-th element is the result of applying the \code{predict} method to the k-th fitted model in \code{object$fits}.
#' @section Author(s):
#' Francis K.C. Hui <francis.hui@anu.edu.au>.
#' @importFrom stats na.pass
#' @examples
#' X <- spider$x
#' abund <- spider$abund
#'
#' # Example 1: Simple example
#' myfamily <- "negative.binomial"
#' # Fit models including all covariates are linear terms, but exclude for bare sand
#' fit0 <- stackedsdm(abund, formula_X = ~. -bare.sand, data = X, family = myfamily, ncores=2) 
#' predict(fit0, type = "response")
#'
#'\donttest{
#' # Example 2: Funkier example where Species are assumed to have different distributions
#' abund[,1:3] <- (abund[,1:3]>0)*1 # First three columns for presence absence
#' myfamily <- c(rep(c("binomial"), 3),
#'        rep(c("negative.binomial"), 5),
#'        rep(c("tweedie"), 4)
#'        )
#' fit0 <- stackedsdm(abund, formula_X = ~ bare.sand, data = X, family = myfamily, ncores=2)
#' predict(fit0, type = "response")
#'}
#' @export predict.stackedsdm
#' @export 
predict.stackedsdm <- function(object, newdata = NULL, type = "link", se.fit = FALSE, na.action = na.pass, ...) {
  num_spp <- length(object$family)
  if(is.null(newdata))
    newdata <- object$data
  
  if(length(type) == 1)
    type <- rep(type, num_spp)
  
  out_preds <- vector("list", num_spp)
  names(out_preds) <- colnames(object$y) 
  for(j in 1:num_spp) {
    if(object$family[j] %in% c("gaussian", "poisson", "Gamma", "binomial", "negative.binomial")) 
      type[j] <- match.arg(type[j], choices = c("response", "link"))
    if(object$family[j] %in% c("zipoisson","zinegative.binomial","ztpoisson","ztnegative.binomial")) 
      type[j] <- match.arg(type[j], choices = c("response", "prob", "count", "zero"))
    if(object$family[j] %in% c("ordinal")) 
      type[j] <- match.arg(type[j], choices = c("prob", "class", "cum.prob", "linear.predictor"))
    
    if(object$family[j] %in% c("gaussian", "poisson", "Gamma", "binomial", "negative.binomial","beta","ztpoisson","ztnegative.binomial","zipoisson","zinegative.binomial","ordinal")) {
      out_preds[[j]]  <- predict(object$fits[[j]]$fit, newdata = newdata, type = type[j], se.fit = se.fit, na.action = na.action)
    }
    if(object$family[j] %in% c("exponential")) {
      out_preds[[j]]  <- predict(object$fits[[j]]$fit, newdata = newdata, type = type[j], se.fit = se.fit, na.action = na.action, dispersion = 1)
    }
    if(object$family[j] == "tweedie") {
      out_preds[[j]]  <- predict(object$fits[[j]]$fit, newdata = newdata, type = type[j], se.fit = se.fit, na.action = na.action)
    }
  }
  if(!se.fit){
    out_preds=data.frame(out_preds)
    colnames(out_preds)= colnames(object$y)  
  }
  
  return(out_preds)
}


#' Calculate residuals from a stackedsdm object
#'
#' @param object An object of class \code{stackedsdm};
#' @param type Determined what type of residuals to calculate. The current options include Dunn-Smyth residuals (default; "dunnsmyth"), raw response residuals ("response") or probability integral transform residuals ("PIT");
#' @param seed For Dunn-Smyth and PIT residuals applied to discrete responses, random jittering is added, and the seed can be used to seed to jittering.
#' @param ... not used
#' @section Details:
#' Calculated the residuals from \code{stackedsdm} object. 
#' @return A matrix of residuals
#' @section Author(s):
#' Francis K.C. Hui <francis.hui@anu.edu.au>.
#' @export residuals.stackedsdm
#' @export 
#' @import stats
#' @importFrom tweedie ptweedie  
#' @examples
#' X <- spider$x
#' abund <- spider$abund
#'
#' # Example 1: Simple example
#' myfamily <- "negative.binomial"
#' # Example 1: Funkier example where Species are assumed to have different distributions
#' # Fit models including all covariates are linear terms, but exclude for bare sand
#' fit0 <- stackedsdm(abund, formula_X = ~. -bare.sand, data = X, family = myfamily, ncores=2) 
#' residuals(fit0)
#'
#' # Example 2: Funkier example where Species are assumed to have different distributions
#' abund[,1:3] <- (abund[,1:3]>0)*1 # First three columns for presence absence
#' myfamily <- c(rep(c("binomial"), 3),
#'               rep(c("negative.binomial"), (ncol(abund)-3)))
#' fit0 <- stackedsdm(abund, formula_X = ~ bare.sand, data = X, family = myfamily, ncores=2)
#' residuals(fit0)
residuals.stackedsdm <- function(object, type = "dunnsmyth", seed = NULL, ...) {
  
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
  
  type <- match.arg(type, choices = c("response", "dunnsmyth", "PIT"))
  num_units <- nrow(object$y)
  num_spp <- ncol(object$y)
  
  out <- object$y - object$fitted
  
  if(type == "response")
    out <- out
  
  if(type %in% c("dunnsmyth", "PIT")) {
    
    for(j in 1:ncol(out)) {
      if(object$family[j] %in% c("gaussian")) 
        out[,j] <- pnorm(object$y[,j], mean = object$fitted[,j], sd = sqrt(object$fits[[j]]$params$dispparam))
      if(object$family[j] %in% c("poisson"))
        out[,j] <- runif(length(object$y[,j]), min = ppois(object$y[,j]-1, lambda = object$fitted[,j]), max = ppois(object$y[,j], lambda = object$fitted[,j]))
      if(object$family[j] %in% c("Gamma")) 
        out[,j] <- pgamma(object$y[,j], scale = object$fits[[j]]$params$dispparam * object$fitted[,j], shape = 1/object$fits[[j]]$params$dispparam)
      if(object$family[j] %in% c("exponential")) 
        out[,j] <- pgamma(object$y[,j], scale = object$fitted[,j], shape = 1)
      if(object$family[j] %in% c("binomial")) {
        if(length(object$trial_size) == 1)
          cw_trial_size <- rep(object$trial_size, num_units)
        if(is.matrix(object$trial_size))
          cw_trial_size <- object$trial_size[,j]
        out[,j] <- runif(num_units, min = pbinom(object$y[,j]-1, size = cw_trial_size, prob = object$fitted[,j]), 
                         max = pbinom(object$y[,j], size = cw_trial_size, prob = object$fitted[,j]))
      }
      if(object$family[j] %in% c("negative.binomial")) 
        out[,j] <- runif(num_units, min = pnbinom(object$y[,j]-1, mu = object$fitted[,j], size = 1/object$fits[[j]]$params$dispparam),
                         max = pnbinom(object$y[,j], mu = object$fitted[,j], size = 1/object$fits[[j]]$params$dispparam))
      if(object$family[j] %in% c("tweedie")) {
        a <- b <- ptweedie(object$y[,j], mu = object$fitted[,j], phi = object$fits[[j]]$params$dispparam, power = object$fits[[j]]$params$powerparam)
        a[object$y[,j]==0] <- 0
        out[,j] <- runif(num_units, min = a, max = b)
      }
      if(object$family[j] %in% c("beta")) 
        out[,j] <- pbeta(object$y[,j], shape1 = object$fitted[,j]/object$fits[[j]]$params$dispparam, 
                         shape2 = (1-object$fitted[,j])/object$fits[[j]]$params$dispparam)
      # if(object$family[j] %in% c("ztpoisson"))
      #      out[,j] <- runif(num_units, min = pztpois(object$y[,j]-1, mean = object$fitted), max = ptzpois(object$y[,j], mean = object$fitted))
      # if(object$family[j] %in% c("ztnegative.binomial")) {
      #      pred_untrunccount <- predict(object$fits[[j]]$fit, type = "count")
      #      out[,j] <- runif(num_units, min = pztnbinom(object$y[,j]-1, mu = pred_untrunccount, size = 1/object$fits[[j]]$params$dispparam), 
      #           max = pztnbinom(object$y[,j], mu = pred_untrunccount, size = 1/object$fits[[j]]$params$dispparam))
      #      rm(pred_untrunccount)
      #      }
      # if(object$family[j] %in% c("zipoisson")) {
      #      pred_untrunccount <- predict(object$fits[[j]]$fit, type = "count")
      #      pred_zeroinfl <- predict(object$fits[[j]]$fit, type = "zero")
      #      out[,j] <- runif(num_units, min = pzipois(object$y[,j]-1, lambda = pred_untrunccount, pi = pred_zeroinfl), 
      #           max = pzipois(object$y[,j], lambda = pred_untrunccount, pi = pred_zeroinfl))
      #      rm(pred_untrunccount, pred_zeroinfl)
      #      }
      # if(object$family[j] %in% c("zinegative.binomial")) {
      #      pred_untrunccount <- predict(object$fits[[j]]$fit, type = "count")
      #      pred_zeroinfl <- predict(object$fits[[j]]$fit, type = "zero")
      #      out[,j] <- runif(num_units, min = pzinbinom(object$y[,j]-1, mu = pred_untrunccount, size = 1/object$fits[[j]]$params$dispparam, pi = pred_zeroinfl), 
      #           max = pzinbinom(object$y[,j]-1, mu = pred_untrunccount, size = 1/object$fits[[j]]$params$dispparam, pi = pred_zeroinfl))
      #      rm(pred_untrunccount, pred_zeroinfl)
      #      }
      if(object$family[j] %in% c("ordinal")) {
        pred_cumprobs <- predict(object$fits[[j]]$fit, type = "cum.prob")
        out[,j] <- runif(num_units, min = pred_cumprobs$cprob2, max = pred_cumprobs$cprob1)
        rm(pred_cumprobs)
      }
    }
    
  }
  
  if(type == "dunnsmyth")
    out <- qnorm(out)
  
  return(out)
}


#' Plot residuals of stackedsdm.
#'
#' @param x is a stackedsdm object.
#' @param ... not used
#' @export plot.stackedsdm
#' @export
#' @examples
#' abund <- spider$abund
#' spider_mod <- stackedsdm(abund,~1, data = spider$x) 
#' plot(spider_mod)
plot.stackedsdm = function(x,...) {
  scatter.smooth(residuals(x)~fitted(x),
                 xlab = "Fitted values",
                 ylab = "Dunn-Smyth residuals")
  abline(h=0,col="red")
}
