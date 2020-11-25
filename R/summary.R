#' Summarising the models fitted for conditional moments
#'
#' Summary methods for the class \code{conditional_moments}. Function gives
#' the summary of Generalised Additive Models fitted for conditional means and
#' conditional variances using \code{conditional_moments} function.
#'
#' @param object an object returned from \code{conditional_moments} function
#' @param type type of the moment to return the model summary. This can take one of
#' "mean" or "variance"
#' @param ... ignored
#'
#' @return returns an object of class summary.gam
#'
#' @seealso \code{\link[mgcv]{gam}}, \code{\link{summary}}
#'
#' @export
summary.conditional_moments <- function(object,
                                        type = c("mean", "variance"), ...){
  type <- match.arg(type)

  if(type == "mean"){
    return(summary(object$mean_gam$fit))
  }
  if(type == "variance"){
    return(summary(object$var_gam$fit))
  }
}
