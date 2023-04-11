#' Estimating conditional mean of a time series
#'
#' This function estimates the means  of a time
#' series conditional on a set of  other times series
#' via additive models.
#'
#' @param data a tibble containing all the time series
#' which are uniquely identified by the corresponding
#' Timestamp.
#' @param formula A GAM formula. See \code{\link[mgcv]{formula.gam}}.
#' The details of model specification are given under ‘Details’.
#' @return The function returns an object of class
#' "gam" as described in \code{\link[mgcv]{gamObject}}.
#' @details{ Suppose $x_t$ is a time series where its
#' mean is a function of $z_t$. i.e. $E(x_t|z_t) = m_x(z_t)$.
#' Then $m_x(z_t)$ can be estimated via generalised
#' additive models (GAM). This function uses
#' GAMs implemented in \code{mgcv} package to estimate
#' the conditional means of a time series given a set of
#' time series predictors.}
#'
#' @seealso \code{\link[mgcv]{gam}}
#' @importFrom mgcv gam
#' @importFrom stats as.formula
#' @export
#' @examples
#' data <- NEON_PRIN_5min_cleaned |>
#'   dplyr::filter(site == "upstream") |>
#'   dplyr::select(Timestamp, turbidity, level, conductance, temperature)
#'
#' fit_mean <- data |>
#'   conditional_mean(turbidity ~ s(level, k = 8) +
#'     s(conductance, k = 8) + s(temperature, k = 8))

conditional_mean <- function(data, formula) {
  # vars <- all.vars(formula)
  # z_fac <- data |> Filter(f = is.factor) |> names
  # y <-  vars[1]
  # z_num <- vars[!(vars %in% c(y, z_fac))]

  # if mean knots are null replace with the default in s()
  # if(is.null(knots_mean)){
  #  knots_mean <- rep(-1, length(z_num))
  # }

  # if(rlang::is_empty(z_fac))
  # {
  # formula_new <- substitute(y ~ s(z_num, k = knots_mean ))

  # formula_new <- paste(y, "~", paste("s(", z_num, ", k=", knots_mean,  ")",
  #                                            sep="", collapse="+"),sep="")

  # } else{
  #  formula_new <- paste(y, "~", paste("s(", z_num, ", k=", knots_mean, ")", sep="", collapse="+"),
  #                         "+", paste(z_fac, collapse = " + "),
  #                        sep = " ")
  # }

  mean_gam_fit <- mgcv::gam(
    formula = stats::as.formula(formula),
    data = data
  )

  mean_gam_fit$type <- "conditional_mean"

  class(mean_gam_fit) <- c("conditional_moment", "gam", "glm", "lm")
  return(mean_gam_fit)

  # return(structure(list(mean_gam_fit),
  #                 class = "conditional_moment"))
}
