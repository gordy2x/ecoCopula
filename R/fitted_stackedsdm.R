#' Fitted values from a stackedsdm object
#'
#' @param object An object of class \code{stackedsdm}
#' @param ... Not used
#' @section Details:
#' Extracts the fitted values from \code{stackedsdm} object. 
#' @return A matrix of fitted values.
#' @section Author(s):
#' Francis K.C. Hui <francis.hui@anu.edu.au>.
#' @export 
#' @examples
#' library(mvabund)
#' data(spider)
#' X <- as.data.frame(spider$x)
#' abund <- spider$abund
#'
#' # Example 1: Simple example
#' myfamily <- "negative.binomial"
#' # Example 1: Funkier example where Species are assumed to have different distributions
#' # Fit models including all covariates are linear terms, but exclude for bare sand
#' fit0 <- stackedsdm(abund, formula_X = ~. -bare.sand, data = X, family = myfamily)
#' fitted(fit0)
#'
#' # Example 2: Funkier example where Species are assumed to have different distributions
#' abund[,1:3] <- (abund[,1:3]>0)*1 # First three columns for presence absence
#' myfamily <- c(rep(c("binomial"), 3),
#'               rep(c("negative.binomial"), (ncol(abund)-3)))
#' fit0 <- stackedsdm(abund, formula_X = ~ bare.sand, data = X, family = myfamily)
#' fitted(fit0)
fitted.stackedsdm <- function(object, ...) {
  return(object$fitted)
}