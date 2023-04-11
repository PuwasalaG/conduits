#' @importFrom tibble tibble
NULL

#' Anomaly removed data for water quality variables aggregated at 5-minute intervals
#' from Pringle Creek, Texas.
#'
#' NEON_PRIN_5min_cleaned consists anomaly removed data for water quality variables
#' from upstream and downstream sensors in Pringle Creek in Texas for the period
#' spanning from 2019-07-01 to 2019-12-31 aggregated at 5-minute intervals.
#'
#' @format A data frame with water-quality variables, level and temperature
#' data:
#' \describe{
#' \item{\code{Timestamp}}{Timestamp}
#' \item{\code{site}}{site position}
#' \item{\code{conductance}}{specific conductance}
#' \item{\code{dissolvedOxygen}}{dissolved oxygen}
#' \item{\code{pH}}{pH}
#' \item{\code{chlorophyll}}{chlorophyll}
#' \item{\code{turbidity}}{turbidity}
#' \item{\code{fDOM}}{fDOM}
#' \item{\code{level}}{elevation of surface water}
#' \item{\code{temperature}}{temperature in surface water}
#' }
"NEON_PRIN_5min_cleaned"
