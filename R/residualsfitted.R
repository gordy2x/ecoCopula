# function() {
#      source("auxfnsv0.R")
#      data(spider)
#      data <- spider$x
#      y <- spider$abun
#      y[,1:3] <- (y[,1:3]>0)*1
#      y[,ncol(y)] <- y[,ncol(y)] + 1
#      family <- rep(c("binomial","poisson","negative.binomial","tweedie"), each = 3)
#      family[ncol(y)] <- "ztnegative.binomial"
# 
#      fit0 <- stackedsdm(y, formula_X = ~. -bare.sand, data = data, family = family)
#      }
#  	

#' Fitted values from a stackedsdm object
#'
#' @param object An object of class \code{stackedsdm}
#' @section Details:
#' Extracts the fitted values from \code{stackedsdm} object. 
#' @return A matrix of fitted values.
#' @section Author(s):
#' Francis K.C. Hui <francis.hui@anu.edu.au>.
#' @export 
#' @examples
#' library(tidyverse)
#' data(spider)
#' data <- spider$x %>% 
#'    as.data.frame
#' y <- spider$abun
#'
#' # Example 1: Simple example
#' myfamily <- "negative.binomial"
#' # Example 1: Funkier example where Species are assumed to have different distributions
#'fit0 <- stackedsdm(y, formula_X = ~. -bare.sand, data = data, family = myfamily) # Fit models including all covariates are linear terms, but exclude for bare sand
#' fitted(fit0)
#'
#' # Example 2: Funkier example where Species are assumed to have different distributions
#' y[,1:3] <- (y[,1:3]>0)*1 # First three columns for presence absence
#' y[,ncol(y)] <- y[,ncol(y)] + 1 # Last column for zero truncated NB
#' myfamily <- rep(c("binomial","poisson","negative.binomial","tweedie"), each = 3)
#' myfamily[ncol(y)] <- "ztnegative.binomial"
#' fit0 <- stackedsdm(y, formula_X = ~. -bare.sand, data = data, family = myfamily)
#' fitted(fit0)
fitted.stackedsdm <- function(object) {
     return(object$fitted)
     }


#' Calculate residuals from a stackedsdm object
#'
#' @param object An object of class \code{stackedsdm};
#' @param type Determined what type of residuals to calculate. The current options include Dunn-Smyth residuals (default; "dunnsmyth"), raw response residuals ("response") or probability integral transform residuals ("PIT");
#' @param seed For Dunn-Smyth and PIT residuals applied to discrete responses, random jittering is added, and the seed can be used to seed to jittering.
#' @section Details:
#' Calculated the residuals from \code{stackedsdm} object. 
#' @return A matrix of residuals
#' @section Author(s):
#' Francis K.C. Hui <francis.hui@anu.edu.au>.
#' @export 
#' @examples
#' library(tidyverse)
#' data(spider)
#' data <- spider$x %>% 
#'    as.data.frame
#' y <- spider$abun
#'
#' # Example 1: Simple example
#' myfamily <- "negative.binomial"
#' # Example 1: Funkier example where Species are assumed to have different distributions
#'fit0 <- stackedsdm(y, formula_X = ~. -bare.sand, data = data, family = myfamily) # Fit models including all covariates are linear terms, but exclude for bare sand
#' residuals(fit0)
#'
#' # Example 2: Funkier example where Species are assumed to have different distributions
#' y[,1:3] <- (y[,1:3]>0)*1 # First three columns for presence absence
#' y[,ncol(y)] <- y[,ncol(y)] + 1 # Last column for zero truncated NB
#' myfamily <- rep(c("binomial","poisson","negative.binomial","tweedie"), each = 3)
#' myfamily[ncol(y)] <- "ztnegative.binomial"
#' fit0 <- stackedsdm(y, formula_X = ~. -bare.sand, data = data, family = myfamily)
#' residuals(fit0)
residuals.stackedsdm <- function(object, type = "dunnsmyth", seed = NULL) {
     type <- match.arg(type, choices = c("response", "dunnsmyth", "PIT"))
     num_units <- nrow(object$y)
     num_spp <- ncol(object$y)

     out <- object$y - object$fitted
     
     if(type == "response")
          out <- out

     if(type %in% c("dunnsmyth", "PIT")) {
          set.seed(seed)
          
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
               if(object$family[j] %in% c("ztpoisson"))
                    out[,j] <- runif(num_units, min = pztpois(object$y[,j]-1, mean = object$fitted), max = ptzpois(object$y[,j], mean = object$fitted))
               if(object$family[j] %in% c("ztnegative.binomial")) {
                    pred_untrunccount <- predict(object$fits[[j]]$fit, type = "count")
                    out[,j] <- runif(num_units, min = pztnbinom(object$y[,j]-1, mu = pred_untrunccount, size = 1/object$fits[[j]]$params$dispparam), 
                         max = pztnbinom(object$y[,j], mu = pred_untrunccount, size = 1/object$fits[[j]]$params$dispparam))
                    rm(pred_untrunccount)
                    }
               if(object$family[j] %in% c("zipoisson")) {
                    pred_untrunccount <- predict(object$fits[[j]]$fit, type = "count")
                    pred_zeroinfl <- predict(object$fits[[j]]$fit, type = "zero")
                    out[,j] <- runif(num_units, min = pzipois(object$y[,j]-1, lambda = pred_untrunccount, pi = pred_zeroinfl), 
                         max = pzipois(object$y[,j], lambda = pred_untrunccount, pi = pred_zeroinfl))
                    rm(pred_untrunccount, pred_zeroinfl)
                    }
               if(object$family[j] %in% c("zinegative.binomial")) {
                    pred_untrunccount <- predict(object$fits[[j]]$fit, type = "count")
                    pred_zeroinfl <- predict(object$fits[[j]]$fit, type = "zero")
                    out[,j] <- runif(num_units, min = pzinbinom(object$y[,j]-1, mu = pred_untrunccount, size = 1/object$fits[[j]]$params$dispparam, pi = pred_zeroinfl), 
                         max = pzinbinom(object$y[,j]-1, mu = pred_untrunccount, size = 1/object$fits[[j]]$params$dispparam, pi = pred_zeroinfl))
                    rm(pred_untrunccount, pred_zeroinfl)
                    }
               if(object$family[j] %in% c("ordinal")) {
                    pred_cumprobs <- predict(object$fits[[j]]$fit, type = "cum.prob")
                    out[,j] <- runif(num_units, min = pred_cumprobs$cprob2, max = pred_cumprobs$cprob1)
                    rm(pred_cumprobs)
                    }
               }
          
          set.seed(NULL)
          }
          
     if(type == "dunnsmyth")
          out <- qnorm(out)
     
     return(out)
     }

