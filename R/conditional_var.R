#' Estimating conditional variance of a time series
#'
#' This function estimates the variance of a time
#' series conditional on a set of  other times series
#' via additive models.
#'
#' @param data A tibble containing all the time series
#' which are uniquely identified by the corresponding
#' Timestamp.
#' @param formula An object of class "formula": a symbolic description of the model to be fitted.
#' The details of model specification are given under ‘Details’.
#' @param family the family to be used in conditional variance model. Currently
#' this can take either "Gamma" or "lognormal".
#' @param fit_mean A GAM object return from  \code{\link[conduits]{conditional_mean}}
#' @return The function returns an object of class
#' "gam" as described in \code{\link[mgcv]{gamObject}}.
#' @details{ Suppose $x_t$ is a time series where its
#' variance is a function of $z_t$. i.e. $Var(x_t|z_t) = v_x(z_t)$.
#' Then $v_x(z_t)$can be estimated via generalised
#' additive models (GAM). This function uses GAMs implemented
#' in \code{mgcv} package to estimate the conditional variance
#' of a time series given a set of time series predictors.}
#'
#' @seealso \code{\link[mgcv]{gam}} and  \code{\link[splines]{ns}}.
#' @importFrom mgcv predict.gam
#' @importFrom mgcv gam
#' @importFrom stats Gamma update as.formula
#' @export
#' @examples
#' data <- NEON_PRIN_5min_cleaned |>
#'   dplyr::filter(site == "upstream") |>
#'   dplyr::select(Timestamp, turbidity, level, conductance, temperature)
#'
#' fit_mean <- data |>
#'   conditional_mean(turbidity ~ s(level, k = 8) +
#'     s(conductance, k = 8) + s(temperature, k = 8))
#' \dontrun{
#' fit_var <- data |>
#'   conditional_var(
#'     turbidity ~ s(level, k = 7) + s(conductance, k = 7) + s(temperature, k = 7),
#'     family = "Gamma",
#'     fit_mean = fit_mean
#'   )
#' }
conditional_var <- function(data, formula,
                            family = c("Gamma", "lognormal"), fit_mean) {
  family <- match.arg(family)
  vars <- all.vars(formula)
  y_name <- vars[1]
  y <- data[[y_name]]

  if (family == "Gamma") {
    formula_var <- stats::update(formula, Y_Ey2 ~ .)
  }
  if (family == "lognormal") {
    formula_var <- stats::update(formula, log(Y_Ey2) ~ .)
  }

  # Compute conditional means, squared errors from the x_mean_gam
  data <- data |>
    dplyr::mutate(
      E_Y = as.numeric(mgcv::predict.gam(fit_mean,
        newdata = data
      )),
      Y_Ey2 = (y - .data$E_Y)^2
    )

  if (family == "Gamma") {
    var_gam_fit <- mgcv::gam(
      formula = stats::as.formula(formula_var),
      data = data,
      family = stats::Gamma(link = "log")
    )
  }

  if (family == "lognormal") {
    var_gam_fit <- mgcv::gam(
      formula = stats::as.formula(formula_var),
      data = data,
      family = stats::gaussian()
    )
  }

  var_gam_fit$type <- "conditional_var"

  class(var_gam_fit) <- c("conditional_moment", "gam", "glm", "lm")
  return(var_gam_fit)
}
