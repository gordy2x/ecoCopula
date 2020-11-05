# Hidden check and family functions

check_family <- function(family, y) {
     if(length(family) != ncol(y) & length(family) != 1) 
          stop("Number of elements in family must either one or the number of columns in y.") 
     if(length(family) == 1) 
          complete_family <- rep(family, ncol(y))
     if(length(family) == ncol(y)) 
          complete_family <- family
     complete_family <- match.arg(complete_family, choices = c("gaussian", "negative.binomial", "poisson", "binomial", "tweedie", "Gamma",
          "exponential", "beta", "ordinal", "ztpoisson", "ztnegative.binomial", "zipoisson", "zinegative.binomial"), several.ok = TRUE)
     if(length(complete_family) != ncol(y))
          stop("At least one of the elements in family is not supported in current version of stackedsdm!")

     if(any(complete_family == "ordinal")) {
          if(sum(y[, complete_family == "ordinal", drop = FALSE] == 0) > 0) 
               stop("For ordinal data, please shift minimum level to 1.")
          }

     if(any(complete_family %in% c("ztpoisson", "ztnegative.binomial"))) {
          if(sum(y[, complete_family == "ztpoisson", drop=FALSE] < 1) > 0)
                    stop("For zero truncated count data, all values have to be greater than or equal to 1.")
          if(sum(y[, complete_family == "ztnegative.binomial", drop=FALSE] < 1) > 0)
                    stop("For zero truncated count data, all values have to be greater than or equal to 1.")
          }
          
     return(complete_family)
     }

     
check_formula_X <- function(formula_X, data) {
     ## Remove any response if put into formula. Using this bare bones approach as it works if you have dots in the formula
     if(length(as.character(formula_X)) == 3)
          stop("formula_X should not have anything on the left hand side of the tilde sign `~'.")
     
     oldform <- stats::as.formula(stats::terms(formula_X, data = data))
     newform <- stats::update.formula(oldform, resp ~ .)

     return(newform)
     }


## Create a negative binomial family with quadratic mean-variance relationship
nb2 <- function() {
	link <- "log"
	linkfun <- function(mu) 
          return(log(mu))
	linkinv <- function(eta) 
          return(pmax(exp(eta), .Machine$double.eps))
	mu.eta <- function(eta) 
          return(pmax(exp(eta), .Machine$double.eps))
	variance <- function(mu, phi) 
          return(pmax(mu+phi*mu^2, .Machine$double.eps))
  
	structure(list(family = "negative.binomial", linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance, name = link), class = "family")
	}
	
	
## Tweedie distribution with log-link
tweedielogfam <- function() {
     link <- "log"
     linkfun <- function(mu) 
          return(log(mu))
     linkinv <- function(eta) 
          return(pmax(exp(eta), .Machine$double.eps))
     mu.eta <- function(eta) 
          return(pmax(exp(eta), .Machine$double.eps))
     variance <- function(mu, power, phi) 
          return(pmax(phi * mu^power, .Machine$double.eps))

     structure(list(family = "tweedie", linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance, link = link), class = "family")
     }

     
## Truncated negative binomial with log-link
# ztnegative.binomial <- function() {
#      link <- "log"
#      linkinv <- function(eta, phi) {
#           return(mean_ztnbinom(mu = pmax(exp(eta), .Machine$double.eps), size = 1/phi))
#           }
#      variance <- function(mu, phi) 
#           return(var_ztnbinom(mu = mu, size = 1/phi))
# 
#      structure(list(family = "ztnegative.binomial", linkinv = linkinv, variance = variance, link = link), class = "family")
#      }
# 
#   
