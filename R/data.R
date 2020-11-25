#' @importFrom tibble tibble
NULL

#' Raw data for water quality variables from Pringle Creek, Texas.
#'
#' NEON_PRIN_waq consists raw data for water quality variables from upstream
#' and downstream sensors in Pringle Creek in Texas for the period
#' spanning from 2019-07-01 to 2019-12-31 at
#' 1-minute intervals.
#'
#' @source This data are retrieved from NEON Aquatic Instrument Data.
#' For further details, see \url{https://data.neonscience.org/data-products/DP1.20288.001}
#'
#'
#' @format A data frame with six water-quality variables with associated
#' anomaly flags:
#' \describe{
#' \item{\code{Timestamp}}{Timestamp}
#' \item{\code{horizontalPosition}}{site position}
#' \item{\code{conductance}}{specific conductance}
#' \item{\code{dissolvedOxygen}}{dissolved oxygen}
#' \item{\code{pH}}{pH}
#' \item{\code{turbidity}}{turbidity}
#' \item{\code{chlorophyll}}{chlorophyll}
#' \item{\code{fDOM}}{fDOM}
#' \item{\code{conductanceAnomalyFlag}}{anomaly flag for conductance}
#' \item{\code{dissolvedOxygenAnomalyFlag}}{anomaly flag for dissolved oxygen}
#' \item{\code{pHAnomalyFlag}}{anomaly flag for pH}
#' \item{\code{turbidityAnomalyFlag}}{anomaly flag for turbidity}
#' \item{\code{chlorophyllAnomalyFlag}}{anomaly flag for chlorophyll}
#' \item{\code{fDOMAnomalyFlag}}{anomaly flag for fDOM}
#' \item{\code{turbidity_wiperAnomalyFlag}}{wiper anomaly flag for turbidity}
#' }
"NEON_PRIN_waq"


#' Raw data for level and temperature variables from Pringle Creek, Texas.
#'
#' NEON_PRIN_level_temp_5min consists raw data for level and temperature variables
#' from upstream and downstream sensors in Pringle Creek in Texas for the period
#' spanning from 2019-07-01 to 2019-12-31 at
#' 5-minute intervals.
#'
#' @source This data are retrieved from NEON Aquatic Instrument Data.
#' For further details, see \url{https://data.neonscience.org/data-products/DP1.20016.001}
#' and \url{https://data.neonscience.org/data-products/DP1.20053.001}.
#'
#'
#' @format A data frame with the following columns:
#' \describe{
#' \item{\code{Timestamp}}{Timestamp}
#' \item{\code{site}}{site position}
#' \item{\code{level}}{elevation of surface water}
#' \item{\code{temperature}}{temperature in surface water}
#' }
"NEON_PRIN_level_temp_5min"


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
#'}
"NEON_PRIN_5min_cleaned"


