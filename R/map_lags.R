#' Computing upstream variables according to the estimated time lag
#'
#' This function will compute the lag of upstream variables according to the estimated time
#' lag returned from \code{estimate_dt} function.
#'
#' @param data a tibble containing all downstream variables, upstream variables, estimated time lag
#' and corresponding Timestamp
#' @param up bared/unquoted name(s) of upstream variable(s)
#' @param down bared/unquoted name(s) of downstream variable(s)
#' @param lag bared/unquoted name corresponds to the estimated time lag
#'
#' @return a tibble with Timestamp, downstream variables and upstream lagged variables corresponds to
#' each timestamp
#'
#' @author Puwasala Gamakumara
#'
#' @export
#'
map_lags <- function(data, up, down, lag){

  names_up <- names(tidyselect::eval_select(dplyr::enquo(up), data))
  names_down <- names(tidyselect::eval_select(dplyr::enquo(down), data))
  names_lag <- names(tidyselect::eval_select(dplyr::enquo(lag), data))

  data <- data %>%
    dplyr::select(.data$Timestamp, {{up}}, {{down}}, {{lag}})

  L_min <- data %>% dplyr::pull({{lag}}) %>% min(na.rm = T)
  L_max <- data %>% dplyr::pull({{lag}}) %>% max(na.rm = T)

  Time <- data$Timestamp

  Time_leads <- purrr::map(L_min:L_max, ~dplyr::lead(data$Timestamp, n = .x))
  names(Time_leads) <- c(L_min:L_max)

  Time_leads <- Time_leads %>%
    tibble::as_tibble() %>%
    dplyr::mutate(Timestamp = Time)

  df_time_lead <- data %>%
    dplyr::select(.data$Timestamp, {{lag}}) %>%
    dplyr::left_join(Time_leads, by = "Timestamp")

  df_ref_time_lead <- df_time_lead %>%
    tidyr::pivot_longer(cols = -c(.data$Timestamp, {{lag}}),
                  names_to = "lead", values_to = "time_d") %>%
    dplyr::select(.data$Timestamp, .data$lead, .data$time_d) %>%
    dplyr::rename(max_lag = .data$lead) %>%
    dplyr::mutate(max_lag = as.double(.data$max_lag))

  # preparing downstream data

  df_downstream <- data %>%
    dplyr::select(.data$Timestamp, {{down}}, {{lag}}) %>%
    dplyr::rename(max_lag = {{lag}})


  # preparing upstream data according to the time delay - we take the average if more than one
  # value correspond to the same time_d.

  up_dt <- paste(names_up, "_dt", sep = "")

  df_upstream <- data %>%
    dplyr::select(.data$Timestamp, {{up}}, {{lag}}) %>%
    dplyr::rename(max_lag = {{lag}}) %>%
    dplyr::inner_join(df_ref_time_lead, by = c("Timestamp", "max_lag")) %>%
    dplyr::group_by(.data$time_d) %>%
    dplyr::summarise_at(dplyr::vars({{up}}), ~mean(.)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$time_d)

  names(df_upstream) <- c("time_d", up_dt)

  # joining upstream and downstream data

  data_final <- df_downstream %>%
    dplyr::left_join(df_upstream, by = c("Timestamp" = "time_d")) %>%
    dplyr::select(.data$Timestamp, tidyselect::everything())

  return(data = data_final)

}
