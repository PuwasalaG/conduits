## code to prepare `NEON_PRIN` dataset goes here

# Code adopted from download-NEON-AIS-data.rmd tutorial

library(neonUtilities)
library(ggplot2)
library(dplyr)
library(padr)
library(stringr)

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


## Splitting waq for upstream and downstream

waq_up <- waq_instantaneous %>%
  filter(horizontalPosition == 101)

waq_down <- waq_instantaneous %>%
  filter(horizontalPosition == 102)

EOS_5min_up <- EOS_5_min %>%
  filter(horizontalPosition == 110)


EOS_5min_down <- EOS_5_min %>%
  filter(horizontalPosition == 132)

TSW_5min_up <- TSW_5min %>%
  filter(horizontalPosition == 101)

TSW_5min_down <- TSW_5min %>%
  filter(horizontalPosition == 102)



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
waq_strm_cols <- base::gsub("dissolvedOxygenSat","dissolvedOxygenSaturation",waq_strm_cols)

waq_final_qf_cols <- variables_20288 %>%
  filter(str_detect(fieldName, "FinalQF")) %>%
  pull(fieldName)


##-- Agggregating into 15 mints interval --##


waq_15min_up <- waq_up %>%
  padr::thicken(interval = "15 min", by = "endDateTime", colname = "roundedTime",
                rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(waq_strm_cols), mean, na.rm = TRUE)

waq_finalqf_15min_up <- waq_up %>%
  padr::thicken(interval = "15 min", by = "endDateTime", colname = "roundedTime",
                rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(waq_final_qf_cols), max, na.rm = TRUE)


waq_15min_down <- waq_down %>%
  padr::thicken(interval = "15 min", by = "endDateTime", colname = "roundedTime", rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(waq_strm_cols), mean, na.rm = TRUE)

waq_finalqf_15min_down <- waq_down %>%
  padr::thicken(interval = "15 min", by = "endDateTime", colname = "roundedTime", rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(waq_final_qf_cols), max, na.rm = TRUE)

# In EOS data we choose only the surfacewaterElevMean variable

EOS_15min_up <- EOS_5min_up %>%
  padr::thicken(interval = "15 min", by = "endDateTime", colname = "roundedTime",
                rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(surfacewaterElevMean), mean, na.rm = TRUE)

EOS_15min_down <- EOS_5min_down %>%
  padr::thicken(interval = "15 min", by = "endDateTime", colname = "roundedTime",
                rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(surfacewaterElevMean), mean, na.rm = TRUE)



# Finally we can merge waq data with the EOS data


all_15min_up <- EOS_15min_up %>%
  left_join(waq_15min_up, by = "roundedTime") %>%
  mutate(site = "110")

all_15min_down <- EOS_15min_down %>%
  left_join(waq_15min_down, by = "roundedTime") %>%
  mutate(site = "132")


# Joining upstream and downstream together

NEON_PRIN_15min <- bind_rows(all_15min_up, all_15min_down)

save(NEON_PRIN_15min, file = "data-raw/NEON_PRIN_15min.rda")




##-- Agggregating into 5 mints interval --##

class(waq_up$endDateTime)

waq_5min_up <- waq_up %>%
  padr::thicken(interval = "5 min", by = "endDateTime", colname = "roundedTime", rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(waq_strm_cols), mean, na.rm = TRUE)

waq_finalqf_5min_up <- waq_up %>%
  padr::thicken(interval = "5 min", by = "endDateTime", colname = "roundedTime", rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(waq_final_qf_cols), max, na.rm = TRUE)



# Aggregating downstream data

class(waq_down$endDateTime)

waq_5min_down <- waq_down %>%
  padr::thicken(interval = "5 min", by = "endDateTime", colname = "roundedTime", rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(waq_strm_cols), mean, na.rm = TRUE)

waq_finalqf_5min_down <- waq_down %>%
  padr::thicken(interval = "5 min", by = "endDateTime", colname = "roundedTime", rounding = "down") %>%
  group_by(roundedTime) %>%
  summarise_at(vars(waq_final_qf_cols), max, na.rm = TRUE)


# Finally we can merge waq data with the EOS data for comparisons


all_5min_up <- list(waq_5min_up, EOS_5min_up, TSW_5min_up) %>%
  purrr::reduce(left_join, by = c("roundedTime" = "endDateTime")) %>%
  mutate(site = "up")

all_5min_down <- list(waq_5min_down, EOS_5min_down, TSW_5min_down) %>%
  purrr::reduce(left_join, by = c("roundedTime" = "endDateTime")) %>%
  mutate(site = "down")


# Joining upstream and downstream together

NEON_PRIN_5min <- bind_rows(all_5min_up, all_5min_down)

NEON_PRIN_5min <- NEON_PRIN_5min %>%
  select(roundedTime, site, surfacewaterElevMean, specificConductance,
         dissolvedOxygen, dissolvedOxygenSaturation, pH, chlorophyll, turbidity, fDOM,
         surfWaterTempMean)

save(NEON_PRIN_5min, file = "data-raw/NEON_PRIN_5min.rda")



usethis::use_data(NEON_PRIN_15min, overwrite = TRUE)
usethis::use_data(NEON_PRIN_5min, overwrite = TRUE)
