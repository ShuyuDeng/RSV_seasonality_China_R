rm(list=ls())
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(aweek)
library(lubridate)
library(wktmo)
source("descriptiveFun.R")
source("weatherFun.R")

# 0. Default----
n.sample <- 100
lags <- c(0, 7, 14, 21, 28)

# 1. Read data----
load("data/rsvdata/ChinaRSVSeasonality_before2020.RData")  # RSV data
load("data/weatherdata/weather_station_2019.RData")        # Weather station info

# 2. Clean RSV data----
str(metadata)
str(monthlydata)
str(weeklydata)
names(metadata)[names(metadata) == "CaseDefinition"] <- "AllCaseDefinitions"        # avoid using same variable name
names(monthlydata)[names(monthlydata) == "Total_year_report"] <- "N_year"
names(weeklydata)[names(weeklydata) == "Total_year_report"] <- "N_year"

# monthly data
monthlydata$StartDate <- as.Date(monthlydata$StartDate)
monthlydata <- left_join(monthlydata, metadata)
monthlydata2 <- monthlydata %>%
  filter(Group == "Main") %>%                # avoid summing duplicate data
  group_by(SID, Latitude, Longitude) %>%
  slice(if (sum(N_year) > 1) {
          which(N_year != 1)                 # prefer to choose multi-year combined data         
        } else {
          seq_len(n())                       # or choose all one-year data
        }) %>%
  ungroup()
monthlydata3 <- monthlydata2 %>% 
  pivot_longer(cols = starts_with(c("M0", "M1")),
               names_to = "MonthSeq",
               names_prefix = "M",
               values_to = "RSV") %>%
  mutate(Month = month(StartDate + months(as.numeric(MonthSeq) - 1)))
# weekly data
weeklydata <- left_join(weeklydata, metadata)
weeklydata2 <- weeklydata %>% 
  rowwise() %>%                                                                                         # group by row
  mutate(StartDate = dateFromWeek(year = YearFrom, wkIndex = WeekIndex, wkMethod = "ISO")[1]) %>%       # impute NA
  ungroup() %>%
  dplyr::select(-WeekIndex) %>%
  filter(Group == "Main") %>% 
  group_by(SID, Latitude, Longitude) %>%
  slice(if (sum(N_year) > 1) {
          which(N_year != 1)
        } else {
          seq_len(n())
        }) %>% 
  ungroup()
# week to month
weeklydata3 <- weeklydata2 %>%
  mutate(numVector = pmap(dplyr::select(., starts_with("W")), 
                          ~unname(unlist(c(...))))) %>%       # convert to numeric vectors
  mutate(wkToMonth = pmap(list(numVector, StartDate), 
                          ~weekToMonth(wkdata = na.omit(..1), datStart = ..2, wkMethod = "startDat"))) %>%
  unnest(wkToMonth) %>%                   # unnest a list-column of data frames into rows and columns
  filter(between(as.Date(paste0(yearMonth, "-01")), 
                 as.Date(paste0(YearFrom, "-", MonthFrom, "-01")), 
                 as.Date(paste0(YearFrom, "-", MonthFrom, "-01")) + months(11))) %>%
  mutate(Month = month(ym(yearMonth))) %>% 
  mutate(RSV = round(value, digits = 0))
# combined data
RSVdata <- rbind(monthlydata3[, c("SID", "Location",  "Province", "Latitude", "Longitude", "Month", "RSV", "N_year")],
                 weeklydata3[, c("SID", "Location",  "Province", "Latitude", "Longitude", "Month", "RSV", "N_year")])
RSVdata2 <- RSVdata %>%
  group_by(Location, Province, Latitude, Longitude, Month) %>%
  summarise(SumRSV = sum(RSV), SumYear = sum(N_year))

# only choose sites with at least 100 cases for analysis
table_summary_case <- RSVdata2 %>%
  group_by(Location, Province, Latitude, Longitude) %>%
  summarise(TotalRSV = sum(SumRSV))
quantile(table_summary_case$TotalRSV)
RSVdata3 <- table_summary_case %>%
  filter(TotalRSV >= n.sample) %>%
  dplyr::select(-TotalRSV) %>%
  left_join(RSVdata2)
rm(monthlydata2, monthlydata3, weeklydata2, weeklydata3)

# 3. Define RSV season----
RSVdata4 <- do.call(rbind, by(RSVdata3, RSVdata3[,c("Latitude", "Longitude")], getAAPEpi, RSVn = "SumRSV", threshold = 75))
RSVdata4 <- RSVdata4 %>%
  mutate(Site = paste(round(Latitude, 6), round(Longitude, 6), sep = "_")) %>%
  left_join(metadata[c("Latitude", "Longitude", "Location", "Province")] %>% distinct(Latitude, Longitude, .keep_all = TRUE))
table_summary2 <- RSVdata4 %>% 
  group_by(Latitude, Longitude) %>%
  summarise(TotalRSV = sum(SumRSV))

# 4. Province-level RSV activity----
# only retain provinces with at least 100 cases
observe_provincedata <- RSVdata2 %>%
  left_join(metadata[c("Latitude", "Longitude", "Location", "Province")] %>% 
              distinct(Latitude, Longitude, .keep_all = TRUE)) %>%
  group_by(Province, Month) %>%
  summarise(SumRSV = sum(SumRSV, na.rm = TRUE),
            TotalYear = sum(SumYear, na.rm = TRUE)) %>%
  mutate(TotalRSV = sum(SumRSV, na.rm = TRUE)) %>%
  filter(TotalRSV >= n.sample) # remove Tianjin（29）、Ningxia（49例）

unique_province <- observe_provincedata %>%
  distinct(Province, TotalYear, TotalRSV)

# Description of RSV activity by province
observe_province <- do.call(rbind, 
                            by(observe_provincedata, observe_provincedata$Province, getAAPEpi, RSVn = "SumRSV", threshold = 75))

# 4. Extract weather data----
nrow(unique(metadata[c("Latitude", "Longitude")]))
by(metadata, metadata[c("Latitude", "Longitude")], get_weather_all)     # download weather data in all stations within 200km

# 5. Clean weather data----
# Calculate monthly daily data: 01/01-12/31
weatherdata <- do.call(rbind, by(RSVdata4, RSVdata4[c("Latitude", "Longitude")], get_weather_site, redo = FALSE))
# Calculate annual mean-centred values for each variable
weatherdata2 <- left_join(weatherdata, 
                          weatherdata %>% group_by(Latitude, Longitude) %>% 
                          summarise(TEMPa = mean(TEMP, na.rm = TRUE),
                                    DEWPa = mean(DEWP, na.rm = TRUE),
                                    WDSPa = mean(WDSP, na.rm = TRUE),
                                    PRCPa = mean(PRCP, na.rm = TRUE),
                                    RHa = mean(RH, na.rm = TRUE),
                                    AHa = mean(AH, na.rm = TRUE)))

weatherdata2$TEMPc <- weatherdata2$TEMP - weatherdata2$TEMPa
weatherdata2$DEWPc <- weatherdata2$DEWP - weatherdata2$DEWPa
weatherdata2$WDSPc <- weatherdata2$WDSP - weatherdata2$WDSPa
weatherdata2$PRCPc <- weatherdata2$PRCP - weatherdata2$PRCPa
weatherdata2$RHc <- weatherdata2$RH - weatherdata2$RHa
weatherdata2$AHc <- weatherdata2$AH - weatherdata2$AHa

# Calculate monthly mean data with different lags
weatherdata3 <- lapply(lags, get_weather_lag, dataset = weatherdata2)
names(weatherdata3) <- sprintf("Lag_%02d", lags)
rm(weatherdata, weatherdata2)

# 6. Merge RSV data and weather data----
RSVdata4$Month <- as.factor(RSVdata4$Month)
dataFinal <- purrr::map(weatherdata3, ~ left_join(RSVdata4, .))

# 7. 2019 China weather data----
load("data/weatherdata/weather_China/weather_China_2019.RData")
load("data/weatherdata/weather_China/weather_Macao_2019.RData")
weather_China_2019 <- rbind(weather_China_2019, weather_Macao_2019)
rm(weather_Macao_2019)

# Add province for each weather station
weather_China_2019 <- weather_China_2019 %>%
  left_join(weatherStation2019 %>% 
              dplyr::select(!c(Index, LONGITUDE, LATITUDE)), by = "STNID")

# 8. Clean RSV data for predict year-on-year variations----
# monthlydata
monthlydata_2 <- monthlydata %>%
  arrange(SID, SiteID, YearFrom) %>%
  filter(N_year == 1 & GroupID == 1) %>%
  group_by(SID, SiteID) %>%
  mutate(Continuous = YearFrom == lag(YearFrom) + 1 | is.na(lag(YearFrom)),
         BlockID = cumsum(!Continuous)) %>%
  group_by(SID, SiteID, BlockID) %>%
  mutate(BlockLength = n()) %>%
  ungroup() %>%
  filter(BlockLength >= 3) %>%
  dplyr::select(-Continuous, -BlockID, -BlockLength)  

monthlydata_3 <- monthlydata_2 %>% 
  pivot_longer(cols = starts_with(c("M0", "M1")),
               names_to = "MonthSeq",
               names_prefix = "M",
               values_to = "RSV") %>%
  mutate(Month = month(StartDate + months(as.numeric(MonthSeq) - 1)),
         Year = if_else(Month < MonthFrom, YearTo, YearFrom))

# weeklydata
weeklydata_2 <- weeklydata %>%
  rowwise() %>%                                                                                         # group by row
  mutate(StartDate = dateFromWeek(year = YearFrom, wkIndex = WeekIndex, wkMethod = "ISO")[1]) %>%       # impute NA
  ungroup() %>%
  arrange(SID, SiteID, YearFrom) %>%
  filter(N_year == 1 & GroupID == 1) %>%
  group_by(SID, SiteID) %>%
  mutate(Continuous = YearFrom == lag(YearFrom) + 1 | is.na(lag(YearFrom)),
         BlockID = cumsum(!Continuous)) %>%
  group_by(SID, SiteID, BlockID) %>%
  mutate(BlockLength = n()) %>%
  ungroup() %>%
  filter(BlockLength >= 3) %>%
  dplyr::select(-Continuous, -BlockID, -BlockLength)  

# week to month
weeklydata_3 <- weeklydata_2 %>%
  mutate(numVector = pmap(dplyr::select(., starts_with("W")), 
                          ~unname(unlist(c(...))))) %>%       # convert to numeric vectors
  mutate(wkToMonth = pmap(list(numVector, StartDate), 
                          ~weekToMonth(wkdata = na.omit(..1), datStart = ..2, wkMethod = "startDat"))) %>%
  unnest(wkToMonth) %>%                   # unnest a list-column of data frames into rows and columns
  filter(between(as.Date(paste0(yearMonth, "-01")), 
                 as.Date(paste0(YearFrom, "-", MonthFrom, "-01")), 
                 as.Date(paste0(YearFrom, "-", MonthFrom, "-01")) + months(11))) %>%
  mutate(Month = month(ym(yearMonth)), 
         Year = if_else(Month < MonthFrom, YearTo, YearFrom)) %>% 
  mutate(RSV = round(value, digits = 0))
# combined data
RSVdata_2 <- rbind(monthlydata_3[, c("SID", "SiteID", "Latitude", "Longitude", "Location", "Province", 
                                     "YearFrom", "MonthFrom", "YearTo", "MonthTo",
                                     "Year", "Month", "RSV", "N_year")],
                   weeklydata_3[, c("SID", "SiteID", "Latitude", "Longitude", "Location", "Province", 
                                    "YearFrom", "MonthFrom", "YearTo", "MonthTo", 
                                    "Year", "Month",  "RSV", "N_year")])

# remove duplicate
RSVdata_2 <- RSVdata_2 %>%
  filter(SID != "L0948")

RSVdata_2$StudySiteID <- paste(RSVdata_2$SID, RSVdata_2$SiteID, sep="_")
length(unique(RSVdata_2$SID)) # 43 studies
length(unique(RSVdata_2$StudySiteID))  # 45 study-sites
