#' Visualising the models fitted for conditional moments
#'
#' A function for visualising the fitted generalised additive models for
#' conditional means and variances from \code{conditional_moments} function
#'
#' @param object an object returned from \code{conditional_moments} function
#' @param type the type of moment want to visualise. This can take one of
#' "mean" or "variance"
#' @param title if TRUE returning the plot with a title
#' @param ... ignored
#'
#' @return returns plots visualising fitted gam model for conditional
#' mean or variance vs each predictor.
#'
#' @importFrom ggplot2 ggplot geom_point aes geom_line geom_ribbon ylab geom_rug
#'
#' @seealso \code{\link[mgcv]{plot.gam}}, \code{\link[visreg]{visreg}}
#'
#' @author Puwasala Gamakumara
#'
#' @export
autoplot.conditional_moments <- function(object,
                                         type = c("mean", "variance"),
                                         title = TRUE,
                                         ...){
  type <- match.arg(type)
  x <- object$x
  z_numeric <- object$z_numeric
  z_factors <- object$z_factors

  mean_fit <- object$mean_gam$fit
  var_fit <- object$var_gam$fit

  names_x <- names(tidyselect::eval_select(dplyr::enquo(x), object$data_conditional_moments))
  names_z_numeric <- names(tidyselect::eval_select(dplyr::enquo(z_numeric),
                                                   object$data_conditional_moments))
  names_z_factors <- names(tidyselect::eval_select(dplyr::enquo(z_factors),
                                                   object$data_conditional_moments))

  p <- length(names_z_numeric) + length(names_z_factors)

  if(p <= 2){
    n_col <- 2
  }

  if(p > 2){
    n_col <- 3
  }


  if(type == "mean"){

    data_vis <- object$mean_gam$data %>%
      dplyr::select(.data$Timestamp, {{x}}, {{z_numeric}}, {{z_factors}}) %>%
      tidyr::drop_na()
    summary_mean <- object$mean_gam$fit %>% summary
    mean_edf <- summary_mean$edf %>% round(digits = 3)

    # getting each component of the linear predictor from the fitted model
    fv_mean <- mgcv::predict.gam(mean_fit, type = "terms")

    if(is.null(z_factors)){
      colnames(fv_mean) <- paste("s_", names_z_numeric, sep = "")
    }else{
      colnames(fv_mean) <- c(paste("s_", names_z_factors, sep = ""),
                             paste("s_", names_z_numeric, sep = ""))
    }

    # residuals
    resid_mean <- mgcv::residuals.gam(mean_fit, type = "working")

    # combining linear predictor and residuals
    resid_mean <- cbind(fv_mean, resid_mean) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(Timestamp = data_vis$Timestamp)

    ##-- visualising fitted mean model
    plot_means <- list()
    data <- list()

    for (i in 1:length(names_z_numeric)) {
      ##-- preparing data for visualisation for numerical variables
      data_vis_z <- data_vis %>%
        dplyr::mutate_at(names_z_numeric[-i],
                  ~median(.,na.rm = T))
      if(!is.null(z_factors)){
        data_vis_z <- data_vis_z %>%
          dplyr::mutate_at(names_z_factors, ~names(which.max(table(.))))
      }

      # predicted values
      predict_fit_mean <- mgcv::predict.gam(mean_fit, newdata = data_vis_z,
                                      se.fit = T) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(LI = as.numeric(.data$fit - .data$se.fit),
               UI = as.numeric(.data$fit + .data$se.fit),
               Timestamp = data_vis_z$Timestamp)

      s_name_zi <- paste("s_", names_z_numeric[i],
                         sep = "")
      data_vis_z <- list(data_vis_z, resid_mean, predict_fit_mean) %>%
        purrr::reduce(dplyr::left_join, by = "Timestamp") %>%
        dplyr::mutate(partial_resid = resid_mean + !!rlang::sym(s_name_zi))

      name_zi <- names_z_numeric[i]
      edf_zi <- mean_edf[i]
      range_zi <- data_vis_z %>% dplyr::pull(!!rlang::sym(name_zi))

      plot_means[[i]] <- data_vis_z %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x = !!rlang::sym(name_zi),
                                         y = .data$partial_resid),
                            color = "#666666", size = 0.5) +
        ggplot2::geom_line(color = "#0099FF",
                           ggplot2::aes(x = !!rlang::sym(name_zi),
                                        y = .data$fit),
                           size = 1) +
        ggplot2::geom_ribbon(ggplot2::aes(x = !!rlang::sym(name_zi),
                                          y = .data$fit,
                                          ymin = .data$LI,
                                          ymax = .data$UI),
                    alpha = 0.5,
                    fill = "grey70") +
        ggplot2::ylab(paste("f(",name_zi, ",", edf_zi, ")", sep = "")) +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()(range_zi))

      data[[i]] <- data_vis_z %>%
        dplyr::select(.data$Timestamp)

    }

    if(title==TRUE){
      gridExtra::grid.arrange(grobs = plot_means,
                                           ncol = n_col,
                                           top = paste("conditional mean of",
                                                       names_x, sep = " "))
    }else{
      gridExtra::grid.arrange(grobs = plot_means,
                              ncol = n_col)
    }

    plot_list <- plot_means
  }


  if(type == "variance"){

    data_vis <- object$var_gam$data %>%
      dplyr::select(.data$Timestamp, {{x}}, .data$X_Ex2, {{z_numeric}}, {{z_factors}}) %>%
      tidyr::drop_na()
    summary_var <- object$var_gam$fit %>% summary
    var_edf <- summary_var$edf %>% round(digits = 3)
    # getting each component of the linear predictor from the fitted model
    fv_var <- mgcv::predict.gam(var_fit, type = "terms")

    if(is.null(z_factors)){
      colnames(fv_var) <- paste("s_", names_z_numeric, sep = "")
    }else{
      colnames(fv_var) <- c(paste("s_", names_z_factors, sep = ""),
                            paste("s_", names_z_numeric, sep = ""))
    }

    # residuals
    resid_var <- mgcv::residuals.gam(var_fit, type = "working")

    # combining linear predictor and residuals
    resid_var <- cbind(fv_var, resid_var) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(Timestamp = data_vis$Timestamp)

    ##-- visualising fitted variance model
    plot_var <- list()

    for (i in 1:length(names_z_numeric)) {
      ##-- preparing data for visualisation for numerical variables
      data_vis_z <- data_vis %>%
        dplyr::mutate_at(names_z_numeric[-i],
                  ~median(.,na.rm = T))
      if(!is.null(z_factors)){
        data_vis_z <- data_vis_z %>%
          dplyr::mutate_at(names_z_factors, ~names(which.max(table(.))))
      }

      # predicted values
      predict_fit_var <- mgcv::predict.gam(var_fit, type = "link",
                                           newdata = data_vis_z,
                                           se.fit = T) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(LI = as.numeric(.data$fit - .data$se.fit),
                      UI = as.numeric(.data$fit + .data$se.fit),
                      Timestamp = data_vis_z$Timestamp)

      s_name_zi <- paste("s_", names_z_numeric[i],
                         sep = "")
      data_vis_z <- list(data_vis_z, resid_var, predict_fit_var) %>%
        purrr::reduce(dplyr::left_join, by = "Timestamp") %>%
        dplyr::mutate(partial_resid = exp(resid_var) + !!rlang::sym(s_name_zi))

      name_zi <- names_z_numeric[i]
      edf_zi <- var_edf[i]
      range_zi <- data_vis_z %>% dplyr::pull(!!rlang::sym(name_zi))

      plot_var[[i]] <- data_vis_z %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x = !!rlang::sym(name_zi),
                                         y = .data$partial_resid),
                            color = "#666666", size = 0.5) +
        ggplot2::geom_line(color = "#0099FF",
                           ggplot2::aes(x = !!rlang::sym(name_zi),
                                        y = (.data$fit)),
                           size = 1) +
        ggplot2::geom_ribbon(ggplot2::aes(x = !!rlang::sym(name_zi),
                                          y = (.data$fit),
                                          ymin = .data$LI,
                                          ymax = .data$UI),
                             alpha = 0.5,
                             fill = "grey70") +
        ggplot2::ylab(paste("f(",name_zi,",",edf_zi, ")", sep = "")) +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()(range_zi))

    }

    if(title==TRUE){
      gridExtra::grid.arrange(grobs = plot_var,
                              ncol = n_col,
                              top = paste("conditional variance of",
                                          names_x, sep = " "))
    }else{
      gridExtra::grid.arrange(grobs = plot_var,
                              ncol = n_col)
    }

    plot_list <- plot_var
  }

  invisible(plot_list)
}


#' Visualising the models fitted for conditional cross-correlations
#'
#' A function for visualising the fitted generalised additive models for
#' conditional means, variances and cross-correlations from
#' \code{conditional_ccf} function
#'
#' @param object an object returned from \code{conditional_ccf} function
#' @param type the type of moment or cross-correlation to visualise.
#' This can take one of "mean", "variance" or "cross-correlation"
#' @param k a vector specifying upto how many lags to be plotted for cross-correlation
#' @param ... further arguments passed to or from other methods
#'
#' @return returns plots visualising fitted gam models for conditional
#' mean or variance vs each predictor.
#'
#' @importFrom ggplot2 ggplot geom_point aes geom_line geom_ribbon ylab geom_rug
#' scale_x_continuous scale_y_continuous
#'
#' @seealso \code{\link[visreg]{visreg}}
#'
#' @author Puwasala Gamakumara
#'
#' @export
autoplot.conditional_ccf <- function(object,
                                     type = c("mean", "variance",
                                              "cross-correlation"),
                                     k = NULL,
                                     ...){

  type <- match.arg(type)
  # x <- object$other$x
  # y <- object$other$y

  z_numeric <- object$other$z_numeric
  z_factors <- object$other$z_factors

  if(is.null(k)){
    k <- object$other$k
  }


  cond_moment_x <- object$data_visualise$conditional_moments$x
  # Only showing the plots for conditional moments of y at lag 1
  cond_moment_y <- object$data_visualise$conditional_moments$y$`k = 1`

  # names_x <- names(eval_select(enquo(x), object$data_ccf))
  names_z_numeric <- names(tidyselect::eval_select(dplyr::enquo(z_numeric), object$data_ccf))
  names_z_factors <- names(tidyselect::eval_select(dplyr::enquo(z_factors), object$data_ccf))

  p <- length(names_z_numeric) + length(names_z_factors)

  if(p <= 2){
    n_col <- 2
  }

  if(p > 2){
    n_col <- 3
  }

  #Conditional mean of x

  oask <- grDevices::devAskNewPage(TRUE)
  on.exit(grDevices::devAskNewPage(oask))


  # visualising conditional means
  if(type == "mean"){
    autoplot.conditional_moments(object = cond_moment_x,
                                 type = "mean")
    autoplot.conditional_moments(object = cond_moment_y,
                                 type = "mean")
  }

  # conditional variance of x

  # visualising conditional variance
  if(type == "variance"){
    autoplot.conditional_moments(object = cond_moment_x,
                                 type = "variance")
    autoplot.conditional_moments(object = cond_moment_y,
                                 type = "variance")
  }

  # visualising conditional cross-correlations

  plot_list <- list()
  plot_list[[1]] <- list()
  plot_list[[2]] <- list()

  if(type == "cross-correlation"){

    for (j in k) {

      ccf_gam_fit_k <- object$data_visualise$conditional_ccf$ccf_gam_fit[[j]]

      data_vis <- ccf_gam_fit_k$data %>%
        tidyr::drop_na()

      # getting each component of the linear predictor from the fitted model
      fv_ccf <- stats::predict.glm(ccf_gam_fit_k, type = "terms")

      if(is.null(z_factors)){
        colnames(fv_ccf) <- paste("s_", names_z_numeric, sep = "")
      }else{
        colnames(fv_ccf) <- c(paste("s_", names_z_factors, sep = ""),
                              paste("s_", names_z_numeric, sep = ""))
      }

      # residuals
      resid_ccf <- stats::residuals.glm(ccf_gam_fit_k, type = "working")

      # combining linear predictor and residuals
      resid_ccf <- cbind(fv_ccf, resid_ccf) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(Timestamp = data_vis$Timestamp)

      plot_link <- list()


      plot_response <- list()

      ## link
      linkfun <- function(mu) {log((1+mu)/(1-mu))}
      ## inverse link
      linkinv <- function(eta) {(exp(eta) - 1)/(exp(eta) + 1)}

      for (i in 1:length(names_z_numeric)) {
        ##-- preparing data for visualisation for numerical variables
        data_vis_z <- data_vis %>%
          dplyr::mutate_at(names_z_numeric[-i],
                    ~median(.,na.rm = T))
        if(!is.null(z_factors)){
          data_vis_z <- data_vis_z %>%
            dplyr::mutate_at(names_z_factors, ~names(which.max(table(.))))
        }

        # predicted values
        predict_fit_ccf_link <- stats::predict.glm(ccf_gam_fit_k, type = "link",
                                            newdata = data_vis_z,
                                            se.fit = T) %>%
          tibble::as_tibble() %>%
          dplyr::rename("fit_link" = .data$fit,
                 "se.fit_link" = .data$se.fit) %>%
          dplyr::mutate(LI_link = as.numeric((.data$fit_link - .data$se.fit_link)),
                 UI_link = as.numeric((.data$fit_link + .data$se.fit_link)),
                 Timestamp = data_vis_z$Timestamp) %>%
          dplyr::select(-.data$residual.scale)

        predict_fit_ccf_res <- stats::predict.glm(ccf_gam_fit_k,
                                           type = "response",
                                           newdata = data_vis_z,
                                           se.fit = T) %>%
          tidyr::as_tibble() %>%
          dplyr::rename("fit_res" = .data$fit,
                 "se.fit_res" = .data$se.fit) %>%
          dplyr::mutate(LI_res = as.numeric((.data$fit_res - .data$se.fit_res)),
                 UI_res = as.numeric((.data$fit_res + .data$se.fit_res)),
                 Timestamp = data_vis_z$Timestamp) %>%
          dplyr::select(-.data$residual.scale)

        s_name_zi <- paste("s_", names_z_numeric[i],
                           sep = "")
        data_vis_z <- list(data_vis_z, resid_ccf, predict_fit_ccf_link,
                           predict_fit_ccf_res) %>%
          purrr::reduce(dplyr::left_join, by = "Timestamp")

        name_zi <- names_z_numeric[i]
        range_zi <- data_vis_z %>% dplyr::pull(!!rlang::sym(name_zi))

        plot_link[[i]] <- data_vis_z %>%
          ggplot() +
          geom_point(aes(x = !!rlang::sym(name_zi), y = .data$XY),
                     color = "#666666", size = 0.5) +
          geom_line(color = "#0099FF", aes(x = !!rlang::sym(name_zi),
                                           y = .data$fit_link),
                    size = 1) +
          geom_ribbon(aes(x = !!rlang::sym(name_zi), y = .data$fit_link,
                          ymin = .data$LI_link, ymax = .data$UI_link),
                      alpha = 0.5,
                      fill = "grey70") +
          scale_x_continuous(breaks = scales::pretty_breaks()(range_zi)) +
          ylab(paste("f(",name_zi,")", sep = ""))

        plot_response[[i]] <- ggplot() +
          geom_line(data = data_vis_z,
                    color = "#0099FF", aes(x = !!rlang::sym(name_zi),
                                           y = .data$fit_res),
                    size = 1) +
          geom_ribbon(data = data_vis_z,
                      aes(x = !!rlang::sym(name_zi), y = .data$fit_res,
                          ymin = .data$LI_res, ymax = .data$UI_res),
                      alpha = 0.5,
                      fill = "grey70") +
          geom_rug(data = data_vis_z %>%
                     dplyr::filter(.data$resid_ccf >= 0),
                   aes(x = !!rlang::sym(name_zi), y = linkinv(.data$resid_ccf)),
                   sides = "t", color = "#666666") +
          geom_rug(data = data_vis_z %>%
                     dplyr::filter(.data$resid_ccf < 0),
                   aes(x = !!rlang::sym(name_zi), y = linkinv(.data$resid_ccf)),
                   sides = "b", color = "#666666") +
          scale_x_continuous(breaks = scales::pretty_breaks()(range_zi)) +
          ylab(bquote(c[.(j)]))

      }

      gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs = plot_link, ncol = n_col,
                               top = "predictor scale"),
                              gridExtra::arrangeGrob(grobs = plot_response, ncol = n_col,
                               top = "response scale"),
                              top = paste("k=", j, sep = ""))
      plot_list[[1]][[j]] <- plot_link
      plot_list[[2]][[j]] <- plot_response
    }

    names(plot_list) <- c("plots_link", "plots_response")

    names(plot_list$plots_link) <- paste("k=", k, sep = "")
    names(plot_list$plots_response) <- paste("k=", k, sep = "")

    invisible(plot_list)


  }

}




#' Visualising the lag time \eqn{d_t}
#'
#' A function to visualise the estimated lag time, \eqn{d_t}, with each
#' of the predictor used is conditional cross-correlation models. The
#' prediction intervals are computed using the Sieve bootstrap method (Buhlmann 1997)
#'
#' @param object an object returned from \code{estimate_dt} function
#' @param interval if TRUE, the bootstrap prediction intervals are computed.
#' This will increase the computation time.Default is setted to FALSE.
#' @param m number of replications when getting boostrap prediction intervals
#' for \eqn{d_t}
#' @param seed setting seed
#' @param ... further arguments passed to or from other methods
#'
#' @return returns plots visualising \eqn{d_t} with each predictor used
#' in conditional cross-correlation functions. Further it returns a list of
#' following components
#' \item{plot_list}{list of plots - ggplot objects}
#' \item{data}{data used to create individual plots}
#'
#'
#' @importFrom ggplot2 ggplot geom_point aes geom_line geom_ribbon ylab
#' geom_rug scale_y_discrete xlab
#'
#' @author Puwasala Gamakumara
#'
#' @references Bühlmann, P. (1997). Sieve bootstrap for time series. Bernoulli, 3(2), 123-148.
#'
#' @export
autoplot.estimate_dt <- function(object,
                                 interval = FALSE,
                                 m = 100,
                                 seed = NULL,
                                 ...){

  # choose the data used to fit conditional-ccf model
  data <- object$conditional_ccf_object$data_visualise$conditional_ccf$data_ccf_fit

  x <- object$conditional_ccf_object$other$x
  y <- object$conditional_ccf_object$other$y
  z_numeric <- object$conditional_ccf_object$other$z_numeric
  z_factors <- object$conditional_ccf_object$other$z_factors
  k <- object$conditional_ccf_object$other$k
  k_min <- object$k_min
  k_max <- object$k_max

  names_x <- names(tidyselect::eval_select(dplyr::enquo(x), object$data_dt))
  names_y <- names(tidyselect::eval_select(dplyr::enquo(y), object$data_dt))
  names_z_numeric <- names(tidyselect::eval_select(dplyr::enquo(z_numeric), data))
  names_z_factors <- names(tidyselect::eval_select(dplyr::enquo(z_factors), data))

  p <- length(names_z_numeric) + length(names_z_factors)

  formula_XY <- object$conditional_ccf_object$formula_gam

  ccf_gam_fit <- object$conditional_ccf_object$data_visualise$conditional_ccf$ccf_gam_fit

  df_correlation = object$conditional_ccf_object$other$df_correlation
  names_XY <- paste("XY", k, "_star", sep = "")

  Time <- data$Timestamp


  # preparing data for each predictor holding other predictors at their medians.
  ### NOTE - This part only consider numerical predictors. Improve this to adupt for factors as well

  list_visualise_data <- list(p)

  # selecting predictors
  data_predictors <- data %>%
    dplyr::select(Timestamp, {{z_numeric}}, {{z_factors}})

  for (i in seq_along(names_z_numeric)) {
    z_i <- data_predictors %>% dplyr::pull(names_z_numeric[i])
    z_i <- seq(from = min(z_i, na.rm = T),
               to = max(z_i, na.rm = T),
               length.out = 300)
    z_i_ <- data %>%
      dplyr::select(Timestamp, dplyr::all_of(names_z_numeric[-i]),
             {{z_factors}}) %>%
      dplyr::mutate_at(names_z_numeric[-i],
                ~ quantile(., probs = 0.5, na.rm = T)) %>%
      dplyr::slice(1:300)

    if(!is.null(z_factors)){
      z_i_ <- z_i_
      dplyr::mutate_at(names_z_factors, ~names(which.max(table(.))))
    }

    z_i_name <- names_z_numeric[i]
    new_data_bs <- tibble::tibble(Timestamp = data$Timestamp[1:300],
                          z_i = z_i)
    colnames(new_data_bs) <- c("Timestamp", z_i_name)

    new_data_bs <- new_data_bs %>%
      dplyr::left_join(z_i_, by = "Timestamp")


    list_visualise_data[[i]] <- list(2)
    list_visualise_data[[i]][[1]] <- new_data_bs

  }

  names(list_visualise_data) <- names_z_numeric

  # computing point estimates of dt
  for (i in 1:p) {

    data_vis_dt <- estimate_dt(object = object$conditional_ccf_object,
                               new_data = list_visualise_data[[i]][[1]])
    data_vis_dt <- data_vis_dt$data_dt %>%
      dplyr::select(.data$Timestamp, names_z_numeric[i], .data$dt,
                    .data$ccf) %>%
      dplyr::rename("max_lag" = .data$dt)

    list_visualise_data[[i]][[2]] <- data_vis_dt
    names(list_visualise_data[[i]]) <- c("new_data", "data_vis_dt")

  }


  ##-- computing interval estimates of dt via Sieve bootstrap method

  if(interval==TRUE){

    ##-- fitting arma model to the residuals from each conditional ccf model

    # predictions from conditional ccf models
    ccf_predictions <- purrr::lift(modelr::gather_predictions, data = data_predictors)(ccf_gam_fit)

    # calculating the residuals from conditional ccf models
    df_ccf_response <- data %>%
      select(all_of(paste0("XY_", k, "_star", sep = "")))
    df_ccf_predictions <- ccf_predictions %>%
      dplyr::mutate(model = factor(.data$model, levels = paste0("k = ", k, sep = ""),
                            labels = paste0("XY", k, "_star_hat", sep = ""))) %>%
      tidyr::spread(key = .data$model, value = .data$pred) %>%
      dplyr::select(-.data$Timestamp, -{{z_numeric}})
    df_ccf_resid <- df_ccf_response - df_ccf_predictions

    # fitting arma models to the gam residuals
    arma_fit <- apply(df_ccf_resid, 2, forecast::auto.arima,
                      approximation = FALSE, stepwise = FALSE,
                      max.q = 0, d = 0, max.order = 20)

    # a function to bootstrap residuals
    func_bootstrap_resid <- function(model, size){

      resid <- stats::residuals(model)
      bs_resid <- sample(resid, size = size, replace = T)

      return(bs_resid)
    }


    N <- nrow(data)
    if(is.null(seed)){
      seed = 1989
    }

    corrl <- corrlink()

    # defining a list for each z
    DF_ccf_bootstrap_z <- list()
    for (i in 1:p) {
      DF_ccf_bootstrap_z[[i]] <- list()
    }

    j = 1
    size_df = 0

    set.seed(seed)
    start_time <- Sys.time()
    while (size_df < m) { # j = 13, 43

      ##-- bootstrapped data by bootstrapping arma residuals
      bs_resid <- purrr::map_df(.x = arma_fit, .f = func_bootstrap_resid,
                         size = N)

      # calculating bootstrap gam residuals
      df_arma_fitted <- purrr::map_df(.x = arma_fit, .f = stats::fitted)
      df_bootstrap_gam_resid <- df_arma_fitted + bs_resid

      # calculating bootstrap data
      bootstrap_data <- df_ccf_predictions + df_bootstrap_gam_resid
      colnames(bootstrap_data) <- paste0("XY", k, "_star", sep = "")

      bootstrap_data <- data_predictors %>%
        dplyr::bind_cols(bootstrap_data)


      ##-- Fitting gam models for conditional ccf at each lag
      bs_ccf_gam_fit <- list()

      for (l in seq_along(names_XY)) {

        data_ccf_gam <- bootstrap_data %>%
          dplyr::select(.data$Timestamp, names_XY[l], {{z_numeric}}) %>%
          dplyr::rename(XY = names_XY[l])

        GLM_out <- tryCatch(stats::glm(formula = stats::as.formula(formula_XY),
                                data = data_ccf_gam,
                                family = stats::gaussian(link = corrl),
                                start = rep(0,(sum(df_correlation))),
                                control = stats::glm.control(maxit = 400)),
                            error = function(e)print("error"))

        suppressWarnings(if(class(GLM_out) == "character") {
          break
        } else {
          bs_ccf_gam_fit[[l]] <- GLM_out
        })



      }

      if(length(bs_ccf_gam_fit) == max(k)){

        names(bs_ccf_gam_fit) <- paste("k = ", min(k):max(k), sep = "")
        for (i in 1:p) {

          bs_df_ccf_max <- purrr::map_dfc(bs_ccf_gam_fit, stats::predict.glm,
                                          newdata = list_visualise_data[[i]][[1]],
                                          type = "response")

          bs_df_ccf_max <- bs_df_ccf_max %>%
            as.data.frame() %>%
            dplyr::mutate(Timestamp = Time[1:300])
          names(bs_df_ccf_max) <- c(k, "Timestamp")

          bs_df_ccf_max <- bs_df_ccf_max %>%
            dplyr::as_tibble() %>%
            tidyr::gather(-.data$Timestamp, key = max_lag,
                          value = ccf) %>%
            dplyr::filter(.data$max_lag %in% seq(k_min, k_max)) %>%
            dplyr::group_by(.data$Timestamp) %>%
            dplyr::slice(which.max(.data$ccf)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(bootstrap_index = rep(j, dplyr::n()))
          # slice(which.max(abs(ccf)))

          # i is for each z and j is for each bootstrap sample
          DF_ccf_bootstrap_z[[i]][[j]] <- bs_df_ccf_max

          # to remove null objects
          DF_ccf_bootstrap_z[[i]] <- purrr::compact(DF_ccf_bootstrap_z[[i]])

        }

      }

      size_df <- length(DF_ccf_bootstrap_z[[1]])
      j = j+1
      # print(c(j, size_df))


    }
    end_time <- Sys.time()

    # calculating bootstrap intervals
    for (i in 1:p) {

      DF_ccf_bootstrap <- dplyr::bind_rows(DF_ccf_bootstrap_z[[i]])

      DF_ccf_bootstrap_CI <- DF_ccf_bootstrap %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(max_lag = as.numeric(.data$max_lag)) %>%
        dplyr::group_by(.data$Timestamp) %>%
        dplyr::summarise(dt_95_LI = round(stats::quantile(.data$max_lag, probs = c(0.025))),
                         dt_95_UI = round(stats::quantile(.data$max_lag, probs = c(0.975))),
                         ccf_95_LI = round(stats::quantile(.data$ccf, probs = c(0.025))),
                         ccf_95_UI = round(stats::quantile(.data$ccf, probs = c(0.975))),
                         dt_80_LI = round(stats::quantile(.data$max_lag, probs = c(0.1))),
                         dt_80_UI = round(stats::quantile(.data$max_lag, probs = c(0.9))),
                         ccf_80_LI = round(stats::quantile(.data$ccf, probs = c(0.1))),
                         ccf_80_UI = round(stats::quantile(.data$ccf, probs = c(0.9))))

      list_visualise_data[[i]][[2]] <- list_visualise_data[[i]][[2]] %>%
        dplyr::left_join(DF_ccf_bootstrap_CI, by = "Timestamp") %>%
        dplyr::mutate(max_lag = factor(.data$max_lag, levels = k_min:k_max))

    }

  }


  plot_list <- list()

  for (i in 1:p) {

    if(interval==TRUE){
      z_i <- names_z_numeric[i]
      plot_list[[i]] <- list_visualise_data[[i]]$data_vis_dt %>%
        ggplot(aes(x = !!rlang::sym(z_i) , y = .data$max_lag)) +
        geom_point() +
        geom_ribbon(aes(ymin = .data$dt_95_LI, ymax = .data$dt_95_UI,
                        group = 1),
                    alpha = 0.3, fill = "gray60") +
        geom_ribbon(aes(ymin = .data$dt_80_LI, ymax = .data$dt_80_UI,
                        group = 1),
                    alpha = 0.3, fill = "gray30") +
        scale_y_discrete(limits = factor(k), breaks = seq(2, max(k), 4)) +
        ylab(expression(italic(d[t]))) + xlab(names_z_numeric[i])
    }

    if(interval==FALSE){
      z_i <- names_z_numeric[i]
      plot_list[[i]] <- list_visualise_data[[i]]$data_vis_dt %>%
        dplyr::filter(!is.na(max_lag)) %>%
        ggplot(aes(x = !!rlang::sym(z_i) , y = .data$max_lag)) +
        geom_point() +
        scale_y_discrete(limits = factor(k), breaks = seq(2, max(k), 4)) +
        ylab(expression(italic(d[t]))) + xlab(names_z_numeric[i])
    }

  }


  if(p <= 2){
    n_col <- 2
  }

  if(p > 2){
    n_col <- 3
  }

  gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs = plot_list, ncol = n_col))
  invisible(list(plot_list = plot_list,
                 data = list_visualise_data))

}

##-- User defined link function
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

