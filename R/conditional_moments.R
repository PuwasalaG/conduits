conditional_moments <- function(data, x, z_numeric, z_factors = NULL, knots_mean = NULL,
                                knots_variance = NULL) {

  # data: A tibble consisting Timestamp, X, Y, and Z1,..,Zp
  # p: number of predictors

  names_x <- names(tidyselect::eval_select(dplyr::enquo(x), data))
  names_z_numeric <- names(tidyselect::eval_select(dplyr::enquo(z_numeric), data))

  if(is.null(knots_mean)){
    knots_mean <- rep(3, length(names_z_numeric))
  }

  if(is.null(knots_variance)){
    knots_variance <- rep(3, length(names_z_numeric))
  }

  if(!is.null(dplyr::enquo(z_factors))){
    names_z_factors <- names(tidyselect::eval_select(dplyr::enquo(z_factors), data))
    formula_x_mean <- paste(names_x, "~", paste("s(", names_z_numeric, ", k=", knots_mean,  ")",
                                           sep = "",  collapse = " + "),
                       "+", paste(names_z_factors, collapse = " + "),
                       sep = " ")
    formula_x_var <- paste("X_Ex2 ~", paste("s(", names_z_numeric, ", k=", knots_variance,  ")",
                                            sep = "", collapse = " + "),
                           "+", paste(names_z_factors, collapse = " + "),
                           sep = " ")

  } else{
    names_z_factors <- NULL
    formula_x_mean <- paste(names_x, "~", paste("s(", names_z_numeric, ", k=", knots_mean,  ")",
                                           sep = "",  collapse = " + "), sep = " ")
    formula_x_var <- paste("X_Ex2 ~", paste("s(", names_z_numeric, ", k=", knots_variance,  ")",
                                            sep = "", collapse = " + "), sep = " ")

    }


  p <- length(c(names_z_numeric, names_z_factors))


  ##-- conditional means --##


  x_mean_gam <- mgcv::gam(formula = as.formula(formula_x_mean), data = data)


  ##-- conditional variance --##

  # first computing conditional means, squared errors from the x_mean_gam and y_mean_gam

  data <- data %>%
    dplyr::mutate(E_X = as.numeric(mgcv::predict.gam(x_mean_gam, newdata = data)),
                  X_Ex2 = ({{x}} - E_X)^2)

  x_var_gam <- mgcv::gam(formula = as.formula(formula_x_var), data = data,
                   family = Gamma(link = "log"))

  data <- data %>%
    dplyr::mutate(Var_X = as.numeric(mgcv::predict.gam(x_var_gam, newdata = data,
                                                       type = "response"))) %>%
    select(-X_Ex2) %>%
    dplyr::rename_with(.fn = ~paste(c("E_", "Var_"), names_x, sep = ""),
                       .cols = c(E_X, Var_X))

  return(structure(list(data_conditional_moments = data,
                        x_mean_gam_fit = x_mean_gam,
                        x_var_gam_fit = x_var_gam,
                        x = dplyr::enexpr(x), z_numeric = dplyr::enexpr(z_numeric),
                        z_factors = dplyr::enexpr(z_factors)),
                   class = "conditional_moments"))

}



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
