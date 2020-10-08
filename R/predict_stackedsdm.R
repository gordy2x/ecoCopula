#' Predictions from a stackedsdm object
#'
#' @param object An object of class \code{stackedsdm}
#' @param newdata Pptionally, a data frame in which to look for variables with which to predict.  If omitted, the covariates from the existing dataset are used.
#' @param type The type of prediction required.  This can be supplied as either a single character string, when is applied to all species, or a vector of character strings of the same length as \code{ncol(object$y)} specifying the type of predictions desired for each species. The exact type of prediction allowed depends precisely on the distribution, but for many there is at least `"link"' which is on the scale of the linear predictors, and ‘"response"’ which is on the scale of the response variable. The values of this argument can be abbreviated.
#' @param se.fit Logical switch indicating if standard errors are required.
#' @param na.action Function determining what should be done with missing values in '"newdata"'. The default is to predict \code{NA}..
#' @section Details:
#'  This function simply applies a for loop, cycling through each fitted model from the \code{stackedsdm} object and then attempting to construct the relevant predictions by applying the relevant \code{predict} method. Please keep in mind no formatting is done to the predictions.
#' @return A list where the k-th element is the result of applying the \code{predict} method to the k-th fitted model in \code{object$fits}.
#' @section Author(s):
#' Francis K.C. Hui <francis.hui@anu.edu.au>.
#' @examples
#' X <- as.data.frame(spider$x)
#' abund <- spider$abund
#'
#' # Example 1: Simple example
#' myfamily <- "negative.binomial"
#' # Fit models including all covariates are linear terms, but exclude for bare sand
#' fit0 <- stackedsdm(abund, formula_X = ~. -bare.sand, data = X, family = myfamily) 
#' predict(fit0, type = "response")
#'
#'\dontrun{
#' # Example 2: Funkier example where Species are assumed to have different distributions
#' abund[,1:3] <- (abund[,1:3]>0)*1 # First three columns for presence absence
#' myfamily <- c(rep(c("binomial"), 3),
#'        rep(c("negative.binomial"), 5),
#'        rep(c("tweedie"), 4)
#'        )
#' fit0 <- stackedsdm(abund, formula_X = ~ bare.sand, data = X, family = myfamily)
#' predict(fit0, type = "response")
#'}
#' @export 
predict.stackedsdm <- function(object, newdata = NULL, type = "link", se.fit = FALSE, na.action = na.pass) {
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
     
     
