#' Estimating conditional moments for time series
#'
#' This function estimates the means and variance of a time series conditional on a set of
#' other times series via additive models.
#'
#'
#' @param data a tibble containing all the time series which are uniquely identified by the
#' corresponding Timestamp.
#' @param x variable for which the conditional moments to be estimated. Should be given as bared/unquoted name
#' @param z_numeric numerical variable(s) use as predictors. Should be given as bared/unquoted names and use c()
#' for multiple variables
#' @param z_factors factor variable(s) use as predictors. Should be given as bared/unquoted names.
#' NULL for empty factors or use c() for multiple variables
#' @param family the family to be used in conditional variance model. Currently
#' this can take either "Gamma" or "lognormal".
#' @param knots_mean  a vector specifying the dimension of the basis in the smooth term fitting for
#' each predictor in the GAM for conditional mean of $x$. Each component of the vector should corresponds to each predictor specified in
#' "z_numeric".
#' @param knots_variance a vector specifying the dimension of the basis in the smooth term fitting for
#' each predictor in the GAM for conditional variance of $x$. Each component of the vector should corresponds to each predictor specified in
#' "z_numeric".
#'
#' @return an object of class "conditional_moments" with the following components
#'  \item{data_conditional_moments}{The original tibble appended with the estimated
#'  conditional means and variance of the time series given in $x$}
#'  \item{x_mean_gam_fit}{fitted GAM model for conditional mean of $x$}
#'  \item{x_var_gam_fit}{fitted GAM model for conditional variance of $x$}
#'
#' @details{ Suppose $x_t$ is a time series where its mean and the variance are functions of $z_t$.
#' i.e. $E(x_t|z_t) = m_x(z_t)$ and $Var(x_t|z_t) = v_x(z_t)$. Then $m_x(z_t)$ and $v_x(z_t)$ can be
#' estimated via generalised additive models (GAM). \code{conditional_moments} uses GAMs implemented
#' in \code{mgcv} package to estimate the conditional means and the variance of a time series given
#' a set of time series predictors.}
#'
#'
#' @author Puwasala Gamakumara
#'
#' @seealso \code{\link[mgcv]{gam}}
#'
#' @export
conditional_moments <- function(data, x, z_numeric,
                                family = c("Gamma", "lognormal"),
                                z_factors = NULL,
                                knots_mean = NULL,
                                knots_variance = NULL) {

  family <- match.arg(family)
  names_x <- names(tidyselect::eval_select(dplyr::enquo(x), data))
  names_z_numeric <- names(tidyselect::eval_select(dplyr::enquo(z_numeric), data))
  names_z_factors <- names(tidyselect::eval_select(dplyr::enquo(z_factors), data))

  # if mean knots are null replace with the default in s()
  if(is.null(knots_mean)){
    knots_mean <- rep(-1, length(names_z_numeric))
  }

  # if variance knots are null replace with the default in s()
  if(is.null(knots_variance)){
    knots_variance <- rep(-1, length(names_z_numeric))
  }

  if(!rlang::is_empty(names_z_factors)){

    formula_x_mean <- paste(names_x, "~", paste("s(", names_z_numeric, ", k=", knots_mean,  ")",
                                                sep = "",  collapse = " + "),
                            "+", paste(names_z_factors, collapse = " + "),
                            sep = " ")
    if(family == "Gamma"){
      formula_x_var <- paste("X_Ex2 ~", paste("s(", names_z_numeric, ", k=", knots_variance,  ")",
                                              sep = "", collapse = " + "),
                             "+", paste(names_z_factors, collapse = " + "),
                             sep = " ")
    }

    if(family == "lognormal"){
      formula_x_var <- paste("log(X_Ex2) ~", paste("s(", names_z_numeric, ", k=", knots_variance,  ")",
                                                   sep = "", collapse = " + "),
                             "+", paste(names_z_factors, collapse = " + "),
                             sep = " ")
    }


  } else{
    formula_x_mean <- paste(names_x, "~", paste("s(", names_z_numeric, ", k=", knots_mean,  ")",
                                                sep = "",  collapse = " + "), sep = " ")
    if(family == "Gamma"){
      formula_x_var <- paste("X_Ex2 ~", paste("s(", names_z_numeric, ", k=", knots_variance,  ")",
                                              sep = "", collapse = " + "), sep = " ")
    }
    if(family == "lognormal"){
      formula_x_var <- paste("log(X_Ex2) ~", paste("s(", names_z_numeric, ", k=", knots_variance,  ")",
                                                   sep = "", collapse = " + "), sep = " ")
    }


  }


  p <- length(names_z_numeric) + length(names_z_factors)


  ##-- conditional means --##


  x_mean_gam <- mgcv::gam(formula = stats::as.formula(formula_x_mean), data = data)


  ##-- conditional variance --##

  # first computing conditional means, squared errors from the x_mean_gam and y_mean_gam

  data <- data %>%
    dplyr::mutate(E_X = as.numeric(mgcv::predict.gam(x_mean_gam, newdata = data)),
                  X_Ex2 = ({{x}} - .data$E_X)^2)

  if(family == "Gamma"){
    x_var_gam <- mgcv::gam(formula = stats::as.formula(formula_x_var),
                           data = data,
                           family = stats::Gamma(link = "log"))
  }

  if(family == "lognormal"){
    x_var_gam <- mgcv::gam(formula = stats::as.formula(formula_x_var),
                           data = data,
                           family = stats::gaussian())
  }

  ##-- Computing Var_x --##

  if(family == "Gamma"){
    data_x_cond_moments <- data %>%
      dplyr::mutate(Var_X = as.numeric(mgcv::predict.gam(x_var_gam, newdata = data,
                                                         type = "response"))) %>%
      dplyr::select(-.data$X_Ex2) %>%
      dplyr::rename_with(.fn = ~paste(c("E_", "Var_"), names_x, sep = ""),
                         .cols = c(.data$E_X, .data$Var_X))

  }

  if(family == "lognormal"){
    data_x_cond_moments <- data %>%
      dplyr::mutate(Var_X = exp(as.numeric(mgcv::predict.gam(x_var_gam, newdata = data)))) %>%
      dplyr::select(-.data$X_Ex2) %>%
      dplyr::rename_with(.fn = ~paste(c("E_", "Var_"), names_x, sep = ""),
                         .cols = c(.data$E_X, .data$Var_X))

  }

  return(structure(list(data_conditional_moments = data_x_cond_moments,
                        mean_gam = list(fit = x_mean_gam,
                                        data = data),
                        var_gam = list(fit = x_var_gam,
                                       data = data),
                        x = dplyr::enexpr(x), z_numeric = dplyr::enexpr(z_numeric),
                        z_factors = dplyr::enexpr(z_factors)),
                   class = "conditional_moments"))

}

# ##-- Conditional moments with only Gamma family for variance
# conditional_moments <- function(data, x, z_numeric, z_factors = NULL, knots_mean = NULL,
#                                 knots_variance = NULL) {
#
#   names_x <- names(tidyselect::eval_select(dplyr::enquo(x), data))
#   names_z_numeric <- names(tidyselect::eval_select(dplyr::enquo(z_numeric), data))
#   names_z_factors <- names(tidyselect::eval_select(dplyr::enquo(z_factors), data))
#
#   if(is.null(knots_mean)){
#     knots_mean <- rep(3, length(names_z_numeric))
#   }
#
#   if(is.null(knots_variance)){
#     knots_variance <- rep(3, length(names_z_numeric))
#   }
#
#   if(!rlang::is_empty(names_z_factors)){
#
#     formula_x_mean <- paste(names_x, "~", paste("s(", names_z_numeric, ", k=", knots_mean,  ")",
#                                            sep = "",  collapse = " + "),
#                        "+", paste(names_z_factors, collapse = " + "),
#                        sep = " ")
#     formula_x_var <- paste("X_Ex2 ~", paste("s(", names_z_numeric, ", k=", knots_variance,  ")",
#                                             sep = "", collapse = " + "),
#                            "+", paste(names_z_factors, collapse = " + "),
#                            sep = " ")
#
#   } else{
#     formula_x_mean <- paste(names_x, "~", paste("s(", names_z_numeric, ", k=", knots_mean,  ")",
#                                            sep = "",  collapse = " + "), sep = " ")
#     formula_x_var <- paste("X_Ex2 ~", paste("s(", names_z_numeric, ", k=", knots_variance,  ")",
#                                             sep = "", collapse = " + "), sep = " ")
#
#     }
#
#
#   p <- length(names_z_numeric) + length(names_z_factors)
#
#
#   ##-- conditional means --##
#
#
#   x_mean_gam <- mgcv::gam(formula = stats::as.formula(formula_x_mean), data = data)
#
#
#   ##-- conditional variance --##
#
#   # first computing conditional means, squared errors from the x_mean_gam and y_mean_gam
#
#   data <- data %>%
#     dplyr::mutate(E_X = as.numeric(mgcv::predict.gam(x_mean_gam, newdata = data)),
#                   X_Ex2 = ({{x}} - .data$E_X)^2)
#
#   x_var_gam <- mgcv::gam(formula = stats::as.formula(formula_x_var),
#                          data = data,
#                          family = stats::Gamma(link = "log"))
#
#   data_x_cond_moments <- data %>%
#     dplyr::mutate(Var_X = as.numeric(mgcv::predict.gam(x_var_gam, newdata = data,
#                                                        type = "response"))) %>%
#     dplyr::select(-.data$X_Ex2) %>%
#     dplyr::rename_with(.fn = ~paste(c("E_", "Var_"), names_x, sep = ""),
#                        .cols = c(.data$E_X, .data$Var_X))
#
#   return(structure(list(data_conditional_moments = data_x_cond_moments,
#                         mean_gam = list(fit = x_mean_gam,
#                                         data = data),
#                         var_gam = list(fit = x_var_gam,
#                                        data = data),
#                         x = dplyr::enexpr(x), z_numeric = dplyr::enexpr(z_numeric),
#                         z_factors = dplyr::enexpr(z_factors)),
#                    class = "conditional_moments"))
#
# }



# autoplot.conditional_moments <- function(object, type = c("mean", "variance"), ...){
#
#   type <- match.arg(type)
#   z_numeric <- object$z_numeric
#   z_factors <- object$z_factors
#
#   x_mean_gam_fit <- object$x_mean_gam_fit
#   x_var_gam_fit <- object$x_var_gam_fit
#   data <- object$data_conditional_moments
#
#   p <- length(enquo(z_numeric)) + length(enquo(z_factors))
#
#   if(p <= 2){
#     n_col <- 2
#   }
#
#   if(p > 2){
#     n_col <- 3
#   }
#
#   # oask <- grDevices::devAskNewPage(TRUE)
#   # on.exit(grDevices::devAskNewPage(oask))
#
#   # visualising conditional means
#   if(type == "mean"){
#     par(mfrow = c(ceiling(p/n_col),n_col))
#     visreg::visreg(x_mean_gam_fit, data = data)
#     mtext("GAM models for conditional mean of x", outer = TRUE, cex = 1)
#
#     # visreg::visreg(x_mean_gam_fit, data = data, plot = FALSE) %>%
#     #   purrr::map(function(x){plot(x, gg = TRUE) + theme_bw()}) %>%
#     #   gridExtra::marrangeGrob(nrow = ceiling(p/n_col), ncol = n_col)
#   }
#
#   # visualising conditional variance
#   if(type == "variance"){
#     par(mfrow = c(ceiling(p/n_col),n_col))
#     visreg::visreg(x_var_gam_fit, data = data, gg = TRUE)
#     mtext("GAM models for conditional variance of x", outer = TRUE, cex = 1)
#   }
#
#
# }
