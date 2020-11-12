## code to prepare `NEON_PRIN_waq`, `NEON_PRIN_level_temp_5min`,
## `NEON_PRIN_15min_cleaned` and `NEON_PRIN_5min_cleaned` datasets
## goes here

# Code adopted from download-NEON-AIS-data.rmd tutorial

library(neonUtilities)
library(ggplot2)
library(dplyr)
library(padr)
library(stringr)
library(lubridate)

waq <- neonUtilities::loadByProduct(dpID="DP1.20288.001", site="PRIN",
                                    startdate="2017-04", enddate="2020-04",
                                    package="expanded",
                                    timeIndex = "5",
                                    token = Sys.getenv("NEON_TOKEN"),
                                    check.size = F)

nsw <- neonUtilities::loadByProduct(dpID="DP1.20033.001", site="PRIN",
                                    startdate="2017-04", enddate="2020-04",
                                    package="expanded",
                                    token = Sys.getenv("NEON_TOKEN"),
                                    check.size = F)

esw <- neonUtilities::loadByProduct(dpID="DP1.20016.001", site="PRIN",
                                    startdate="2017-04", enddate="2020-04",
                                    package="expanded",
                                    token = Sys.getenv("NEON_TOKEN"),
                                    check.size = F)

tsw <- neonUtilities::loadByProduct(dpID="DP1.20053.001", site="CARI",
                                    startdate="2017-04", enddate="2020-04",
                                    package="expanded",
                                    timeIndex = "5",
                                    token = Sys.getenv("NEON_TOKEN"),
                                    check.size = F)

list2env(waq, .GlobalEnv)
list2env(nsw, .GlobalEnv)
list2env(esw, .GlobalEnv)
list2env(tsw, .GlobalEnv)

# There two sensor locations for water-quality measures at PRIN. 101-upstream, 102-downstream
waq_instantaneous %>% distinct(horizontalPosition)

# There is one sensor location for surface water at PRIN. 102-downstream
NSW_15_minute %>% distinct(horizontalPosition)

# There are four sensor locations for surface water at PRIN. 131-upstream, 132-downstream
EOS_30_min %>% distinct(horizontalPosition)
EOS_5_min %>% distinct(horizontalPosition)

TSW_5min %>% distinct(horizontalPosition)


# saving level and temperature data

EOS_5_min <- EOS_5_min %>%
  mutate(Timestamp = ymd_hms(endDateTime)) %>%
  select(Timestamp, surfacewaterElevMean, horizontalPosition) %>%
  filter(Timestamp >= ymd("2019-07-01") & Timestamp < ymd("2020-01-01")) %>%
  rename("site" = horizontalPosition,
         "level" = surfacewaterElevMean)

TSW_5min <- TSW_5min %>%
  mutate(Timestamp = ymd_hms(endDateTime)) %>%
  select(Timestamp, surfWaterTempMean, horizontalPosition) %>%
  filter(Timestamp >= ymd("2019-07-01") & Timestamp < ymd("2020-01-01")) %>%
  rename("site" = horizontalPosition,
         "temperature" = surfWaterTempMean)

NEON_PRIN_level_temp_5min <- TSW_5min %>%
  left_join(EOS_5_min) %>%
  mutate(site = recode(site, "101" = "upstream", "102" = "downstream")) %>%
  select(Timestamp, site, level, temperature)

usethis::use_data(NEON_PRIN_level_temp_5min, overwrite = TRUE)


## Splitting data for upstream and downstream

waq_up <- waq_instantaneous %>%
  mutate(Timestamp = ymd_hms(endDateTime)) %>%
  filter(horizontalPosition == 101)

waq_down <- waq_instantaneous %>%
  mutate(Timestamp = ymd_hms(endDateTime)) %>%
  filter(horizontalPosition == 102)

EOS_5min_up <- EOS_5_min %>%
  filter(site == 110)

EOS_5min_down <- EOS_5_min %>%
  filter(site == 132)

TSW_5min_up <- TSW_5min %>%
  filter(site == 101)

TSW_5min_down <- TSW_5min %>%
  filter(site == 102)



# Filtering data between 2019-07-01 and 2019-12-31

waq_down <- waq_down %>%
  filter(Timestamp >= ymd("2019-07-01") & Timestamp < ymd("2020-01-01"))
waq_up <- waq_up %>%
  filter(Timestamp >= ymd("2019-07-01") & Timestamp < ymd("2020-01-01"))



##-- Flagging anomalous points in waq --##

# - we identified these flags by inspecting all the variables
# - also we flag the anomalies corresponds to the wiper in

waq_down <- waq_down %>%
  as_tibble() %>%
  mutate(Timestamp = ymd_hms(endDateTime),
         site = "102",
         turbidityAnomalyFlag =
           case_when(Timestamp %in% ymd_hms("2019-07-19 19:52:59",
                                            "2019-10-24 07:45:00",
                                            "2019-10-25 05:20:00",
                                            "2019-10-29 23:17:30",
                                            "2019-10-30 05:18:30",
                                            "2019-10-30 05:34:30",
                                            "2019-11-24 14:48:00",
                                            "2019-11-25 13:14:00",
                                            "2019-11-20 17:52:00",
                                            "2019-11-26 14:23:00",
                                            "2019-11-27 13:17:00",
                                            "2019-12-03 22:27:00",
                                            "2019-12-04 16:20:00",
                                            "2019-12-04 16:31:00",
                                            "2019-12-07 13:10:00",
                                            "2019-12-17 22:45:00")
                     ~ 1,
                     turbidityRangeQF == 1 ~ 1,
                     is.na(turbidity) ~ as.numeric(NA),
                     TRUE ~ 0),
         pHAnomalyFlag = case_when(pHRangeQF == 1 ~ 1,
                                   is.na(pH) ~ as.numeric(NA),
                                   TRUE ~ 0),
         conductanceAnomalyFlag = case_when(is.na(specificConductance)
                                            ~ as.numeric(NA),
                                            specificConductanceRangeQF == 1
                                            ~ 1,
                                            TRUE ~ 0),
         dissolvedOxygenAnomalyFlag = case_when(is.na(dissolvedOxygen)
                                                ~ as.numeric(NA),
                                                dissolvedOxygenRangeQF == 1
                                                ~ 1,
                                                TRUE ~ 0),
         chlorophyllAnomalyFlag =
           case_when(chlorophyllRangeQF == 1 ~ 1,
                     is.na(chlorophyll) ~ as.numeric(NA),
                     TRUE ~ 0),
         fDOMAnomalyFlag =
           case_when(fDOMRangeQF == 1 ~ 1,
                     is.na(fDOM) ~ as.numeric(NA),
                     TRUE ~ 0)) %>%
  select(Timestamp, horizontalPosition, specificConductance,
         conductanceAnomalyFlag, dissolvedOxygen, dissolvedOxygenAnomalyFlag,
         pH, pHAnomalyFlag, turbidity, turbidityAnomalyFlag, chlorophyll,
         chlorophyllAnomalyFlag, fDOM, fDOMAnomalyFlag)

waq_up <- waq_up %>%
  as_tibble() %>%
  mutate(Timestamp = ymd_hms(endDateTime),
         site = "101",
         turbidityAnomalyFlag =
           case_when(Timestamp %in% ymd_hms("2019-07-30 16:15:30",
                                            "2019-08-06 05:47:30",
                                            "2019-09-26 19:22:30",
                                            "2019-10-01 19:03:30",
                                            "2019-10-11 12:46:30",
                                            "2019-10-12 15:17:30",
                                            "2019-10-16 20:04:30",
                                            "2019-11-13 01:44:00",
                                            "2019-11-16 12:48:00",
                                            "2019-11-18 13:20:00",
                                            "2019-11-20 11:33:00")
                     ~ 1,
                     turbidityRangeQF == 1 ~ 1,
                     is.na(turbidity) ~ as.numeric(NA),
                     TRUE ~ 0),
         pHAnomalyFlag = case_when(Timestamp >= ymd("2019-08-18")
                                   & Timestamp <= ymd("2019-08-21") ~ 1,
                                   pHRangeQF == 1 ~ 1,
                                   is.na(pH) ~ as.numeric(NA),
                                   TRUE ~ 0),
         conductanceAnomalyFlag = case_when(is.na(specificConductance)
                                            ~ as.numeric(NA),
                                            specificConductanceRangeQF == 1
                                            ~ 1,
                                            TRUE ~ 0),
         dissolvedOxygenAnomalyFlag = case_when(is.na(dissolvedOxygen)
                                                ~ as.numeric(NA),
                                                dissolvedOxygenRangeQF == 1
                                                ~ 1,
                                                TRUE ~ 0),
         chlorophyllAnomalyFlag =
           case_when(Timestamp %in% ymd_hms("2019-09-27 08:44:30",
                                            "2019-10-01 22:57:30",
                                            "2019-10-23 07:45:00",
                                            "2019-12-05 09:26:00",
                                            "2019-12-17 16:31:30",
                                            "2019-12-22 04:53:30")
                     ~ 1,
                     chlorophyllRangeQF == 1 ~ 1,
                     is.na(chlorophyll) ~ as.numeric(NA),
                     TRUE ~ 0),
         fDOMAnomalyFlag = as.numeric(NA)) %>%
  select(Timestamp, horizontalPosition, specificConductance,
         conductanceAnomalyFlag, dissolvedOxygen, dissolvedOxygenAnomalyFlag,
         pH, pHAnomalyFlag, turbidity, turbidityAnomalyFlag, chlorophyll,
         chlorophyllAnomalyFlag, fDOM, fDOMAnomalyFlag)


# Wiper flags in turbidity sensor

waq_up <- waq_up %>%
  mutate(turbidity_max = zoo::rollmax(turbidity, k = 5, fill = NA,
                                      align = "right"),
         turbidity_wiperAnomalyFlag = if_else(turbidity == turbidity_max, 1, 0),
         turbidity_wiperAnomalyFlag = if_else(is.na(turbidity_wiperAnomalyFlag), 0,
                                              turbidity_wiperAnomalyFlag),
         turbidityAnomalyFlag = turbidityAnomalyFlag
         + turbidity_wiperAnomalyFlag) %>%
  select(-turbidity_max)



waq_down <- waq_down %>%
  mutate(turbidity_max = zoo::rollmax(turbidity, k = 5, fill = NA,
                                      align = "right"),
         turbidity_wiperAnomalyFlag = if_else(turbidity == turbidity_max, 1, 0),
         turbidity_wiperAnomalyFlag = if_else(is.na(turbidity_wiperAnomalyFlag),
                                              0, turbidity_wiperAnomalyFlag),
         turbidityAnomalyFlag = turbidityAnomalyFlag + turbidity_wiperAnomalyFlag) %>%
  select(-turbidity_max)


NEON_PRIN_waq <- bind_rows(waq_up, waq_down) %>%
  rename("conductance" = specificConductance)

usethis::use_data(NEON_PRIN_waq, overwrite = TRUE)





###-- Data aggregation --###

# We will aggregate data into 15 minutes and 5 minutes intervals.
# EOS has 5 minutes interval data which can be used to aggregate into 15 mint intervals.
# NSW is only available in 15 mint intervals. Therefore the variables in NSW (mostly NO3 data)
# are not available in 5mint invervals.
# We also want to aggregate waq data product also into 15-minutes interval.

# At PRIN, the elevation of surface water sensor is co-located with the water
# quality sonde at horizontalPosition = "101" and "102", meaning the upstream and downstream
# sensor set.

# First we find the water-quality variable names

waq_strm_betaqf_cols <- variables_20288 %>%
  filter(str_detect(fieldName, "BetaQF")) %>%
  pull(fieldName)

waq_strm_cols <- base::gsub("BetaQF","",waq_strm_betaqf_cols)
waq_strm_cols <- purrr::map(c("Sat", "Depth"), str_subset,
                            string = waq_strm_cols, negate = T) %>%
  purrr::reduce(intersect)

waq_final_qf_cols <- variables_20288 %>%
  filter(str_detect(fieldName, "FinalQF")) %>%
  pull(fieldName)
waq_final_qf_cols <- purrr::map(c("Sat", "Depth"), str_subset,
                                string = waq_final_qf_cols,
                                negate = T) %>%
  purrr::reduce(intersect)

# Now remove the anomalies identified in the previos section

waq_up_cleaned <- waq_up %>%
  mutate(turbidity = if_else(turbidityAnomalyFlag == 1, as.numeric(NA),
                             turbidity),
         specificConductance = if_else(conductanceAnomalyFlag == 1,
                                       as.numeric(NA), specificConductance),
         pH = if_else(pHAnomalyFlag == 1, as.numeric(NA),
                      pH),
         dissolvedOxygen = if_else(dissolvedOxygenAnomalyFlag == 1,
                                   as.numeric(NA), dissolvedOxygen),
         chlorophyll = if_else(chlorophyllAnomalyFlag == 1,
                               as.numeric(NA), chlorophyll),
         fDOM = if_else(fDOMAnomalyFlag == 1, as.numeric(NA),
                        fDOM))

waq_down_cleaned <- waq_down %>%
  mutate(turbidity = if_else(turbidityAnomalyFlag == 1, as.numeric(NA),
                             turbidity),
         specificConductance = if_else(conductanceAnomalyFlag == 1,
                                       as.numeric(NA), specificConductance),
         pH = if_else(pHAnomalyFlag == 1, as.numeric(NA),
                      pH),
         dissolvedOxygen = if_else(dissolvedOxygenAnomalyFlag == 1,
                                   as.numeric(NA), dissolvedOxygen),
         chlorophyll = if_else(chlorophyllAnomalyFlag == 1,
                               as.numeric(NA), chlorophyll),
         fDOM = if_else(fDOMAnomalyFlag == 1, as.numeric(NA),
                        fDOM))



##-- Agggregating into 15 mints interval --##


waq_15min_up <- waq_up_cleaned %>%
  padr::thicken(interval = "15 min", by = "Timestamp",
                colname = "roundedTime",
                rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(waq_strm_cols), mean, na.rm = TRUE)

# waq_finalqf_15min_up <- waq_up %>%
#   padr::thicken(interval = "15 min", by = "Timestamp", colname = "roundedTime",
#                 rounding = "down") %>%
#   group_by(roundedTime) %>%
#   summarise_at(vars(waq_final_qf_cols), max, na.rm = TRUE)


waq_15min_down <- waq_down_cleaned %>%
  padr::thicken(interval = "15 min", by = "Timestamp",
                colname = "roundedTime", rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(waq_strm_cols), mean, na.rm = TRUE)

# waq_finalqf_15min_down <- waq_down %>%
#   padr::thicken(interval = "15 min", by = "endDateTime", colname = "roundedTime", rounding = "down") %>%
#   group_by(roundedTime) %>%
#   summarise_at(vars(waq_final_qf_cols), max, na.rm = TRUE)

# In EOS data we choose only the surfacewaterElevMean variable

EOS_15min_up <- EOS_5min_up %>%
  padr::thicken(interval = "15 min", by = "Timestamp", colname = "roundedTime",
                rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(level), mean, na.rm = TRUE)

EOS_15min_down <- EOS_5min_down %>%
  padr::thicken(interval = "15 min", by = "Timestamp", colname = "roundedTime",
                rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(level), mean, na.rm = TRUE)

# surface temperature data
TSW_15min_up <- TSW_5min_up %>%
  padr::thicken(interval = "15 min", by = "Timestamp", colname = "roundedTime",
                rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(temperature), mean, na.rm = TRUE)

TSW_15min_down <- TSW_5min_down %>%
  padr::thicken(interval = "15 min", by = "Timestamp", colname = "roundedTime",
                rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(temperature), mean, na.rm = TRUE)


# Finally we can merge waq data with the EOS data


all_15min_up <- list(waq_15min_up, EOS_15min_up, TSW_15min_up) %>%
  purrr::reduce(left_join, by = "roundedTime") %>%
  mutate(site = "upstream")

all_15min_down <- list(waq_15min_down, EOS_15min_down, TSW_15min_down) %>%
  purrr::reduce(left_join, by = "roundedTime") %>%
  mutate(site = "downstream")


# Joining upstream and downstream together

NEON_PRIN_15min_cleaned <- bind_rows(all_15min_up, all_15min_down) %>%
  rename("conductance" = specificConductance,
         "Timestamp" = "roundedTime") %>%
  select(Timestamp, site, everything())


usethis::use_data(NEON_PRIN_15min_cleaned, overwrite = TRUE)




##-- Agggregating into 5 mints interval --##

class(waq_up$endDateTime)

waq_5min_up <- waq_up_cleaned %>%
  padr::thicken(interval = "5 min", by = "Timestamp", colname = "roundedTime", rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(waq_strm_cols), mean, na.rm = TRUE)

# waq_finalqf_5min_up <- waq_up_cleaned %>%
#   padr::thicken(interval = "5 min", by = "endDateTime", colname = "roundedTime", rounding = "down") %>%
#   group_by(roundedTime) %>%
#   summarise_at(vars(waq_final_qf_cols), max, na.rm = TRUE)



# Aggregating downstream data

class(waq_down$endDateTime)

waq_5min_down <- waq_down_cleaned %>%
  padr::thicken(interval = "5 min", by = "Timestamp", colname = "roundedTime", rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(waq_strm_cols), mean, na.rm = TRUE)

# waq_finalqf_5min_down <- waq_down_cleaned %>%
#   padr::thicken(interval = "5 min", by = "endDateTime", colname = "roundedTime", rounding = "down") %>%
#   group_by(roundedTime) %>%
#   summarise_at(vars(waq_final_qf_cols), max, na.rm = TRUE)


# Finally we can merge waq data with the EOS data for comparisons


all_5min_up <- list(waq_5min_up, EOS_5min_up %>% select(-site),
                    TSW_5min_up %>% select(-site)) %>%
  purrr::reduce(left_join, by = c("roundedTime" = "Timestamp")) %>%
  mutate(site = "upstream")

all_5min_down <- list(waq_5min_down, EOS_5min_down %>% select(-site),
                      TSW_5min_down %>% select(-site)) %>%
  purrr::reduce(left_join, by = c("roundedTime" = "Timestamp")) %>%
  mutate(site = "downstream")


# Joining upstream and downstream together

NEON_PRIN_5min_cleaned <- bind_rows(all_5min_up, all_5min_down) %>%
  rename("conductance" = specificConductance,
         "Timestamp" = "roundedTime") %>%
  select(Timestamp, site, everything())


usethis::use_data(NEON_PRIN_5min_cleaned, overwrite = TRUE)
