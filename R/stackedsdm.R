#' Stacked species regression models, possibly fitted in parallel 
#'
#' @param y A matrix of species responses
#' @param formula_X An object of class \code{formula} representing the relationship to the covariates to be fitted. There should be nothing to the left hand side of the "~" sign.
#' @param data Data frame of the covariates
#' @param family Either a single character vector, in which case all responses are assumed to be from this family, or a vector of character strings of the same length as the number of columns of y. Families as strings and not actual \code{family} class objects. This could be changed though if desired in the future e.g., for custom link functions. Currently, the following families are supported (hopefully properly!): "gaussian", "negative.binomial" (with quadratic mean-variance relationship), "poisson", "binomial" (with logit link), "tweedie", "Gamma" (with log link), "exponential", "beta" (with logit link), "ordinal" (cumulative logit model), "ztpoisson", "ztnegative.binomial", "zipoisson", "zinegative.binomial". 
#' @param trial_size The trial size if any of the responses are binomial. Is either a single number or a matrix with the same dimension as y. If the latter, then all columns that do not correspond to binomial responses are ignored.
#' @param do_parallel Do the separate species model fits in parallel? Defaults to \code{TRUE}
#' @param ncores The number of cores to use if separate the species model fits are done in parallel. If \code{do_parallel = TRUE}, then it defaults to \code{detectCores() - 2}
#' @param trace Print information. This is not actually used currently
#' @section Details:
#' \code{stackedsdm} behaves somewhat like the \code{manyglm} or \code{manyany} function in the package \code{\link{mvabund}}, in the sense that it fits a separate regression to each species response i.e., column of \code{y}. The main difference is that different families can be permitted for each species, which thus allows for mixed responses types. 
#' @return A object of class \code{stackedsdm} with the following components:
#' \code{call} The function call;
#' \code{fits} A list where the j-th element corresponds to the to the fitted model for species j i.e., the j-th column in \code{y};
#; \code{y, formula_X, data, family, trial_size} As per the arguments;
#' \code{linear_predictor} A matrix of the fitted linear predictors
#' \code{fitted} A matrix of the fitted values
#' @section Author(s):
#' Francis K.C. Hui <francis.hui@anu.edu.au>.
#' @import glm2
#' @import mgcv
#' @import mvabund
#' @importFrom  betareg betareg
#' @import ordinal 
#' @import compiler
#' @import doParallel 
#' @import foreach
#' @import stats
#' @export 
#' @examples
#' data(spider)
#' X <- spider$x
#' abund <- spider$abund
#'
#' # Example 1: Simple example
#' myfamily <- "negative.binomial"
#' # Example 1: Funkier example where Species are assumed to have different distributions
#' # Fit models including all covariates are linear terms, but exclude for bare sand
#' fit0 <- stackedsdm(abund, formula_X = ~. -bare.sand, data = X, family = myfamily, ncores = 2) 
#'
#' # Example 2: Funkier example where Species are assumed to have different distributions
#' abund[,1:3] <- (abund[,1:3]>0)*1 # First three columns for presence absence
#' myfamily <- c(rep(c("binomial"), 3),
#'               rep(c("negative.binomial"), (ncol(abund)-3)))
#' fit0 <- stackedsdm(abund, formula_X = ~ bare.sand, data = X, family = myfamily, ncores = 2)
stackedsdm <- function(y, formula_X= ~1, data=NULL, family="negative.binomial", 
                       trial_size = 1, do_parallel = FALSE, 
                       ncores = NULL, trace = FALSE) {
        y <- as.matrix(y)
        if(is.null(colnames(y)))
                colnames(y) <- paste0("resp", 1:ncol(y))
        if(is.null(rownames(y)))
                rownames(y) <- paste0("units", 1:nrow(y))
        num_units <- nrow(y)
        num_spp <- ncol(y)
        
        formula_X <- check_formula_X(formula_X, data = data)
        
        family <- check_family(family = family, y = y)
        if(is.matrix(trial_size)) {
                if(nrow(trial_size) != num_units | ncol(trial_size) != num_spp)
                        stop("trial_size should either be a scalar, or a matrix of the same as y. If the latter, then columns which are not classed as having binomial responses are ignored.")
        }
        
        if(do_parallel) {
                if(is.null(ncores))
                        registerDoParallel(cores = parallel::detectCores()-2)
                if(!is.null(ncores))
                        registerDoParallel(cores = ncores)
        }
        
        
        ##---------------
        ## Set up function to do SDM per species
        ##---------------
        respfit_fn <- function(j) {
                out_params <- list()
                
                if(family[j] %in% c("gaussian", "poisson")) {
                        fit_init <- glm2(formula_X, family = family[j], data = data.frame(resp = y[,j], data))
                        out_params$coefficients <- fit_init$coefficients
                        if(family[j] == "gaussian")
                                out_params$dispparam <- summary(fit_init)$sigma^2
                }
                if(family[j] %in% c("Gamma", "exponential")) {
                        fit_init <- glm2(formula_X, family = Gamma(link = "log"), data = data.frame(resp = y[,j], data))
                        out_params$coefficients <- fit_init$coefficients
                        if(family[j] == "exponential") {
                                out_params$coefficients <- summary(fit_init, disperson = 1)$coefficients[,1]
                        }
                }
                if(family[j] %in% c("binomial")) {
                        if(length(trial_size) == 1)
                                cw_trial_size <- rep(trial_size, num_units)
                        if(is.matrix(trial_size))
                                cw_trial_size <- trial_size[,j]
                        formula_X <- update.formula(formula_X, cbind(resp, trial_size - resp) ~ .)
                        fit_init <- glm2(formula_X, family = "binomial", data = data.frame(resp = y[,j], trial_size = cw_trial_size, data))
                        out_params$coefficients <- fit_init$coefficients
                }
                if(family[j] %in% c("negative.binomial")) {
                        fit_init <- manyglm(formula_X, data = data.frame(resp = y[,j], data), family = "negative.binomial")
                        out_params$coefficients <- fit_init$coefficients
                        out_params$dispparam<- fit_init$phi
                }
                if(family[j] == "tweedie") {
                        fit_init <- gam(formula_X, data = data.frame(resp = y[,j], data), family = mgcv::tw(), method = "ML")
                        out_params$coefficients <- fit_init$coefficients
                        out_params$dispparam <- summary(fit_init)$dispersion
                        out_params$powerparam <- as.numeric(strsplit(strsplit(fit_init$family$family, "p=")[[1]][2], ")")[[1]])
                }
                if(family[j] == "beta") {
                        fit_init <- betareg(formula_X, data = data.frame(resp = y[,j], data), link = "logit") 
                        out_params$coefficients <- fit_init$coefficients$mean
                        out_params$dispparam <- 1/fit_init$coefficients$precision
                }
                # if(family[j] == "ztpoisson") {
                #      fit_init <- countreg::zerotrunc(formula_X, data = data.frame(resp = y[,j], data), dist = "poisson")
                #      out_params$coefficients <- fit_init$coefficients
                #      }
                # if(family[j] == "ztnegative.binomial") {
                #      fit_init <- countreg::zerotrunc(formula_X, data = data.frame(resp = y[,j], data), dist = "negbin")
                #      out_params$coefficients <- fit_init$coefficients
                #      out_params$dispparam <- 1/fit_init$theta
                #      }
                # if(family[j] == "zipoisson") {
                #      fit_init <- countreg::zeroinfl(formula_X, data = data.frame(resp = y[,j], data), dist = "poisson", link = "logit")
                #      out_params$coefficients <- fit_init$coefficients$count
                #      out_params$ziintercept <- fit_init$coefficients$zero
                #      }
                # if(family[j] == "zinegative.binomial") {
                #      fit_init <- countreg::zeroinfl(formula_X, data = data.frame(resp = y[,j], data), dist = "negbin", link = "logit")
                #      out_params$coefficients <- fit_init$coefficients$count
                #      out_params$ziintercept <- fit_init$coefficients$zero
                #      out_params$dispparam <- 1/fit_init$theta
                #      }
                if(family[j] == "ordinal") {
                        fit_init <- clm(formula_X, data = data.frame(resp = ordered(y[,j]), data), link = "logit") 
                        out_params$coefficients <- fit_init$beta
                        out_params$cutoffs <- fit_init$alpha
                }
                
                return(list(params = out_params, fit = fit_init))
        }
        
        respfit_cmpfn <- compiler::cmpfun(respfit_fn)
        rm(respfit_fn)
        
        j=NULL #need this to pass check()
        if(do_parallel)
                all_fits <- foreach(j = 1:num_spp) %dopar% respfit_cmpfn(j = j)          
        if(!do_parallel)
                all_fits <- foreach(j = 1:num_spp) %do% respfit_cmpfn(j = j)          
        
        
        ##-----------------
        ## Format output
        ##-----------------
        out_allfits <- list(call = match.call(), fits = all_fits, y = y, formula_X = formula_X, data = data, family = family, trial_size = trial_size)
        out_allfits$linear_predictor <- sapply(all_fits, function(x) x$fit$linear)
        out_allfits$fitted <- sapply(all_fits, function(x) fitted(x$fit))
        dimnames(out_allfits$fitted)=dimnames(y)
        class(out_allfits) <- "stackedsdm"
        
        rm(all_fits)
        gc()
        return(out_allfits)
}
     
     
#' @export 
mgcv::ldTweedie