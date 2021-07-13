#' Computing conditional cross-correlations at given lags
#'
#' This function computes cross correlation between $x_t$ and $y_{t+k}$ at $k = 1,2,...$
#' conditional on a set of time series $z_t$
#'
#' @param data a tibble containing all the time series which are uniquely identified by the
#' corresponding Timestamp
#' @param x bared/unquoted name of the variable to be considered as $x_t$
#' @param y bared/unquoted name of the variable to be considered as $y_t$
#' @param z_numeric numerical variable(s) use as predictors. Should be given as bared/unquoted names
#' and use c() for multiple variables
#' @param z_factors factor variable(s) use as predictors. Should be given as bared/unquoted names.
#' NULL for empty factors or use c() for multiple variables
#' @param family the family to be used in conditional variance model. Currently
#' this can take either "Gamma" or "lognormal".
#' @param k a vector of lag values at which the cross-correlation needs to be computed. Default is
#' $1:9$
#' @param knots_mean a named list of vectors specifiying the dimension of the basis of the smooth term fitting for
#' each predictor in the conditional mean models for $x_t$ and $y_t$ (see \code{conditional_moments}). The vectors should be named as
#' $x$ and $y$. The components of each vector should correspond to each predictor specified in
#' "z_numeric". By default function fits a 3 dimentional thin plate regression spline for each predictor.
#' @param knots_variance a named list of vectors specifiying the dimension of the basis of the smooth term fitting for
#' each predictor in the conditional variance models for $x_t$ and $y_t$ (see \code{conditional_moments}). The vectors should be named as
#' $x$ and $y$. The components of each vector should correspond to each predictor specified in
#' "z_numeric". By default the function fits a 3 dimentional thin plate regression spline for each predictor.

#' @param df_correlation a vector specifying the degrees of freedom to be considered for each numerical
#' predictor when fitting additive models for conditional cross-correlations. Each component of the
#' vector should corresponds to each predictor specified in "z_numeric". By default the function
#' will fit a natural cubic spline with $3$ degrees of freedom. see \code{\link[splines]{ns}}
#' consider $3$ degrees for each predictor.
#'
#' @return an object of class "conditional_ccf" with the following components
#'  \item{data_ccf}{The original tibble passed to the function appended with the estimated
#'  conditional cross correlation at lags $k = 1,2,...$}
#'  \item{data_visualise}{a list containing the fitted models for conditional means, variance and
#'  cross-correlation at each lag. These data can be used for visualisation and diagnostics of each fitted model}
#'  \item{formula_gam}{formula passed to the \code{\link[stats]{glm}} to fit additive models for
#'  conditional cross-correlation}
#'
#' @details{ Suppose $x_t$ and $y_t$ are conditionally normalised with respect
#' to $z_t$ using \code{conditional_moments}. Then we can estimate the conditional
#' cross-correlation between $x_t$ and $y_t$ at lag $k$, i.e. $r_k = E(x_ty_{t+k}|z_t)$
#' via generalised additive models (GAM). \code{conditional_ccf} uses natural splines implemented
#' in \code{splines} package to estimate the conditional cross-correlations between two
#' time series given a set of time series predictors. Users need not to
#' normalise $x_t$ and $y_t$. The function \code{conditional_ccf} itself will
#' normalise $x_t$ and $y_t$ using \code{conditional_moments}}
#'
#' @author Puwasala Gamakumara
#'
#' @seealso \code{\link[stats]{glm}}, \code{\link[splines]{ns}}
#'
#' @export
conditional_ccf <- function(data, x, y, z_numeric, z_factors,
                            family = c("Gamma", "lognormal"),
                            k = 1:9,
                            knots_mean = NULL, knots_variance = NULL, df_correlation = NULL){

  family <- match.arg(family)
  names_x <- names(tidyselect::eval_select(dplyr::enquo(x), data))
  names_y <- names(tidyselect::eval_select(dplyr::enquo(y), data))
  names_z_numeric <- names(tidyselect::eval_select(dplyr::enquo(z_numeric), data))
  names_z_factors <- names(tidyselect::eval_select(dplyr::enquo(z_factors), data))

  p_numeric <- length(names_z_numeric)
  p <- length(names_z_numeric) + length(names_z_factors)

  if(is.null(knots_mean)){
    knots_mean = list(x = rep(3, p_numeric), y = rep(3, p_numeric))
  }

  if(is.null(knots_variance)){
    knots_variance = list(x = rep(3, p_numeric), y = rep(3, p_numeric))
  }

  if(is.null(df_correlation)){
    df_correlation = rep(3, p_numeric)
  }


  data_NEW <- data %>%
    dplyr::select(.data$Timestamp, {{x}}, {{y}}, {{z_numeric}}, {{z_factors}})


  ##-- Computing conditional moments --##

  data_x_cond_moments <- data_NEW %>%
    dplyr::select(.data$Timestamp, {{x}}, {{z_numeric}}, {{z_factors}})
  data_y_cond_moments <- data_NEW %>%
    dplyr::select(.data$Timestamp, {{y}}, {{z_numeric}}, {{z_factors}})

  #computing conditional moments for x

  cond_moments_x <- conditional_moments(data = data_x_cond_moments,
                                        x = {{x}},
                                        family = family,
                                        z_numeric = {{z_numeric}},
                                        knots_mean = knots_mean$x,
                                        knots_variance = knots_variance$x)

  data_x_cond_moments <- cond_moments_x$data_conditional_moments %>%
    dplyr::select(-{{z_numeric}}, -!!rlang::enexpr(names_z_factors))

  #computing conditional moments for y_t+k at each k

  data_y_k_cond_moments <- list()
  cond_moments_y_out <- list()

  for (i in k) {
    data_y <- data_y_cond_moments %>%
      dplyr::mutate_at(vars({{y}}), dplyr::lead, n=i)

    cond_moments_y_out[[i]] <- conditional_moments(data = data_y,
                                                   x = {{y}},
                                                   family = family,
                                                   z_numeric = {{z_numeric}},
                                                   knots_mean = knots_mean$y,
                                                   knots_variance = knots_variance$y)

    data_y_k_cond_moments[[i]] <- cond_moments_y_out[[i]]$data_conditional_moments %>%
      dplyr::select(-{{z_numeric}}, -!!rlang::enexpr(names_z_factors))


  }


  #computing the conditionally normalised X_t and Y_t+k, i.e., X* and Y*_t+k

  E_X <- paste("E_", names_x, sep = "")
  E_Y <- paste("E_", names_y, sep = "")
  Var_X <- paste("Var_", names_x, sep = "")
  Var_Y <- paste("Var_", names_y, sep = "")

  #for x
  data_x_cond_moments <- data_x_cond_moments %>%
    dplyr::mutate(X_star = as.numeric(({{x}} - !!rlang::sym(E_X))/sqrt(!!rlang::sym(Var_X))))

  #for y_t+k
  for (i in k) {
    data_y_k_cond_moments[[i]] <- data_y_k_cond_moments[[i]] %>%
      dplyr::mutate(Y_star = as.numeric(({{y}} - !!rlang::sym(E_Y))/sqrt(!!rlang::sym(Var_Y))))

  }

  ##-- Computing x*y*_t+k --##

  DF_XY_star <- data_x_cond_moments %>%
    select(Timestamp)

  for (i in k) {
    X_star <- data_x_cond_moments %>% pull(X_star)
    Y_star <- data_y_k_cond_moments[[i]] %>% pull(Y_star)
    DF_XY_star <- DF_XY_star %>%
      dplyr::mutate("XY_{i}_star" := X_star*Y_star)

  }


  data_estim_r_gam <- data_NEW %>%
    dplyr::select(.data$Timestamp, {{z_numeric}},
                  !!rlang::enexpr(names_z_factors)) %>%
    dplyr::left_join(DF_XY_star, by = "Timestamp")


  ##-- Computing cross-correlation at lags k --##

  if(!rlang::is_empty(names_z_factors)){
    formula_XY <- paste("XY ~ - 1 +", paste("splines::ns(", names_z_numeric, ", df=",
                                            df_correlation,  ")", sep = "", collapse = " + "),
                        "+", paste(names_z_factors, collapse = " + "),
                        sep = " ")
  }else{
    z_factors <- NULL
    formula_XY <- paste("XY ~ - 1 +", paste("splines::ns(", names_z_numeric, ", df=",
                                            df_correlation,  ")", sep = "",
                                            collapse = " + "), sep = " ")
  }


  ccf_gam_fit <- list()
  DF_ccf_max <- matrix(0, ncol = length(k), nrow = nrow(data_NEW))
  corrl <- corrlink()
  names_XY <- paste("XY_", k, "_star", sep = "")

  for (i in k) {

    data_ccf_gam <- data_estim_r_gam %>%
      dplyr::select(.data$Timestamp, names_XY[i], {{z_numeric}}, !!rlang::enexpr(names_z_factors)) %>%
      dplyr::rename(XY = names_XY[i])


    ccf_gam_fit[[i]] <- stats::glm(formula = stats::as.formula(formula_XY),
                                   data = data_ccf_gam,
                                   family = stats::gaussian(link = corrl),
                                   start = rep(0,(sum(df_correlation))),
                                   control = stats::glm.control(maxit = 400))

    DF_ccf_max[,i] <- stats::predict.glm(ccf_gam_fit[[i]], newdata = data_NEW, type = "response")

  }

  names(ccf_gam_fit) <- paste("k = ", k, sep = "")
  colnames(DF_ccf_max) <- paste("k = ", k, sep = "")

  DF_ccf_max <- DF_ccf_max %>%
    tibble::as_tibble() %>%
    dplyr::mutate(Timestamp = data_NEW$Timestamp) %>%
    dplyr::select(.data$Timestamp, tidyselect::everything())

  DF_ccf_max <- data %>%
    dplyr::left_join(DF_ccf_max, by = "Timestamp")

  #naming cond_moments_y_out
  names(cond_moments_y_out) <- paste("k = ", k, sep = "")

  return(structure(list(data_ccf = DF_ccf_max,
                        data_visualise = list(conditional_moments = list(x = cond_moments_x,
                                                                         y = cond_moments_y_out),
                                              conditional_ccf = list(ccf_gam_fit = ccf_gam_fit,
                                                                     data_ccf_fit = data_estim_r_gam)),
                        formula_gam = formula_XY,
                        other = list(x = rlang::enexpr(x), y = rlang::enexpr(y),
                                     z_numeric = rlang::enexpr(z_numeric),
                                     z_factors = rlang::enexpr(z_factors),
                                     k = k,
                                     knots_mean = knots_mean,
                                     knots_variance = knots_variance,
                                     df_correlation = df_correlation)),
                   class = "conditional_ccf"))

}

compute_XY_star <- function(data, k){

  Time <- data$Timestamp
  Y_leads <- purrr::map(k, ~dplyr::lead(data$Y_star, n = .x))
  names(Y_leads) <- paste("Y_", k, sep = "")

  Y_leads <- Y_leads %>%
    tibble::as_tibble()

  # Computing X_start*Y_star_leads

  X_star <- data %>% dplyr::pull(X_star)

  func_XY <- function(y, X_star){
    y*X_star
  }

  Y_leads_new <- Y_leads %>%
    dplyr::mutate_all(~ func_XY(., X_star = X_star)) %>%
    dplyr::mutate(Timestamp = Time)

  names_XY <- paste("XY", k, "_star", sep = "")
  colnames(Y_leads_new) <- c(names_XY, "Timestamp")

  # Joining X*Y* and filtering for the complete cases

  data_estim_r_gam <- list(data, Y_leads_new) %>%
    purrr::reduce(dplyr::left_join, by = "Timestamp")

  return(data_estim_r_gam)

}

corrlink <- function() {
  ## link
  linkfun <- function(mu) {log((1+mu)/(1-mu))}
  ## inverse link
  linkinv <- function(eta) {(exp(eta) - 1)/(exp(eta) + 1)}
  ## derivative of invlink wrt eta
  mu.eta <- function(eta) { 2*exp(eta)/(exp(eta) + 1)^2 }
  valideta <- function(eta) TRUE
  link <- "corrlink"
  structure(list(linkfun = linkfun,
                 linkinv = linkinv,
                 mu.eta = mu.eta,
                 valideta = valideta,
                 name = link),
            class = "link-glm")
}
