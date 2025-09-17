library(GSODR)
library(geosphere)
library(imputeTS)

# Download all meteorological data from weather stations within 200 km of each study site across all study years
get_weather_all <- function(site) {
  Latitude <- round(site$Latitude[1], 6)
  Longitude <- round(site$Longitude[1], 6)
  minYearA <- min(site$YearA)
  maxYearB <- max(site$YearB)
  filename <- paste(Latitude, "_", Longitude, ".Rdata", sep = "")
  existing_files <- list.files("data/weatherdata/weather_all/", pattern = "\\.Rdata$", full.names = FALSE)
  # Check if the data is already downloaded
  if (filename %in% existing_files) {
    message("Downloaded: ", Latitude, "_", Longitude)
    return(NULL)
  } else {
    # Extract all weather data <200 km
    station_ids <- nearest_stations(Latitude, Longitude, distance = 200)
    message(Latitude, "_", Longitude)
    weather_all <- suppressWarnings({
      do.call(rbind, lapply(station_ids, function(station_id) {
        message("StationID: ", station_id, " (", which(station_id == station_ids), "/", length(station_ids), ")")
        get_GSOD(minYearA: maxYearB, station = station_id)
        }))
    })
    save(weather_all, file = paste("weather_all/", filename, sep = ""))
    # Remove data and release memory
    rm(weather_all)
    gc()
    return(NULL)
    }
}

# Calculate absolute humidity
calcAH <- function(RH, TEMP) {
  AH <- 6.112*exp((17.67*TEMP)/(TEMP + 243.5))*RH*2.1674/(273.15 + TEMP)
  return(AH)
}

# Clean meteorological data for each study site; if multiple years exist, compute daily monthly averages
get_weather_site <- function(site, redo = FALSE) {
  Latitude <- site$Latitude[1]
  Longitude <- site$Longitude[1]
  message(Latitude, "_", Longitude)
  if(redo) {
  load(file = paste("data/weatherdata/weather_all/", round(Latitude, 6), "_", round(Longitude, 6), ".Rdata", sep = ""))
  weather_site <- weather_all[, c("STNID", "NAME", "LATITUDE", "LONGITUDE", 
                                  "YEARMODA", "TEMP", "DEWP", "WDSP", "PRCP", "PRCP_ATTRIBUTES", "RH")]
  rm(weather_all)
  ## Only include 24h daily Rainfall data; otherwise marked as NA
  weather_site$PRCP[is.na(weather_site$PRCP)] <- 0
  weather_site$PRCP[weather_site$PRCP_ATTRIBUTES!="G" & 
                      weather_site$PRCP_ATTRIBUTES!="F" & 
                      weather_site$PRCP_ATTRIBUTES!="D"] <- NA
  ## Add AH
  weather_site$AH <- calcAH(RH = weather_site$RH, TEMP = weather_site$TEMP)
  ## Only include the nearest weather station for each date
  coordinates <- data.frame(
    YEARMODA = weather_site$YEARMODA,
    Longitude = Longitude,
    Latitude = Latitude,
    LONGITUDE = weather_site$LONGITUDE,
    LATITUDE = weather_site$LATITUDE)
  nearest_stations <- coordinates %>%
    group_by(YEARMODA) %>%
    mutate(Distance = distVincentyEllipsoid(cbind(LONGITUDE, LATITUDE), cbind(Longitude, Latitude))) %>%
    filter(Distance == min(Distance))
  ## Get mean weather data on the same date for multiple years
  weather_site2 <- nearest_stations %>% 
    left_join(weather_site, by = c("YEARMODA", "LONGITUDE", "LATITUDE")) %>% 
    mutate(MODA = format(YEARMODA, "%m-%d")) %>% 
    group_by(MODA) %>%
    summarise(TEMP = mean(TEMP, na.rm = TRUE),
              DEWP = mean(DEWP, na.rm = TRUE),
              WDSP = mean(WDSP, na.rm = TRUE),
              PRCP = mean(PRCP, na.rm = TRUE),
              RH = mean(RH, na.rm = TRUE),
              AH = mean(AH, na.rm = TRUE),
              LATITUDE = mean(LATITUDE, na.rm = TRUE),
              LONGITUDE = mean(LONGITUDE, na.rm = TRUE),
              Latitude = mean(Latitude),
              Longitude = mean(Longitude)) %>%
    filter(MODA != "02-29") %>% 
    mutate(DaySeq = order(MODA))
  write.csv(weather_site2, file = paste("data/weatherdata/weather_site/", round(Latitude, 6), "_", round(Longitude, 6), ".csv", sep = ""), row.names = FALSE)
  return(weather_site2)
  }else{
  weather_site <- read.csv(file = paste("data/weatherdata/weather_site/", round(Latitude, 6), "_", round(Longitude, 6), ".csv", sep = ""))
  return(weather_site)
  }
}

# Obtain lagged meteorological data
get_weather_lag <- function(dataset, lag) {
  weather_lag <- dataset %>%
    mutate(Month = cut(ifelse(DaySeq - lag < 1, DaySeq - lag + 365, DaySeq - lag),
                       breaks = c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365),
                       labels = c(1:12), right = TRUE)) %>%
    group_by(Latitude, Longitude, Month) %>%
    summarise(TEMP = mean(TEMP, na.rm = TRUE),
              DEWP = mean(DEWP, na.rm = TRUE),
              WDSP = mean(WDSP, na.rm = TRUE),
              PRCP = mean(PRCP, na.rm = TRUE),
              RH = mean(RH, na.rm = TRUE),
              AH = mean(AH, na.rm = TRUE),
              TEMPc = mean(TEMPc, na.rm = TRUE),
              DEWPc = mean(DEWPc, na.rm = TRUE),
              WDSPc = mean(WDSPc, na.rm = TRUE),
              PRCPc = mean(PRCPc, na.rm = TRUE),
              RHc = mean(RHc, na.rm = TRUE),
              AHc = mean(AHc, na.rm = TRUE))
  return(weather_lag)
}

# Download all meteorological data for China in 2019 from GSOD
get_weather_China_2019 <- function() {

  # China 2019 weather data: no Macao
  weather_China_2019 <- get_GSOD(years = 2019, country = "China")
  unique(weather_China_2019$Province)
  save(weather_China_2019, file = "data/weatherdata/weather_China/weather_China_2019.RData")

  # Add Macao 2019 weather data
  weather_Macao_2019 <- get_GSOD(years = 2019, country = "MO")
  save(weather_Macao_2019, file = "data/weatherdata/weather_China/weather_Macao_2019.RData")
  rm(weather_China_2019, weather_Macao_2019)
}


# Clean daily meteorological data for each province in 2019
get_weather_province_2019 <- function(weatherdata, redo = FALSE) {
  message(weatherdata$Province[1], weatherdata$name[1])

  if (redo) {

    # Calculate absolute humidity, group by province, compute daily averages in 2019
    weather_province_2019 <- weatherdata %>%
      mutate(AH = calcAH(RH = RH, TEMP = TEMP)) %>%
      group_by(YEARMODA) %>%
      summarise(Province = Province[1],
                DaySeq = YDAY[1],
                TEMP = mean(TEMP, na.rm = TRUE), 
                AH = mean(AH, na.rm = TRUE)) %>%
      mutate(DaySeq = order(YEARMODA))
    missing.percent <- (365 - n_distinct(weather_province_2019$DaySeq))/365 * 100
    
    # If missing data >5%, discard the province's 2019 weather data
    if (missing.percent > 5) {
      return(NULL)
    } else {
      
      # Otherwise, create a complete date sequence and add province info
      completeDates <- seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "day")
      df <- data.frame(Province = weatherdata$Province[1],
                       YEARMODA = completeDates, 
                       DaySeq = as.numeric(format(completeDates, "%j")))
      weather_province_2019 <- df %>%
        left_join(weather_province_2019, by = c("Province", "YEARMODA", "DaySeq"))

      # Interpolate missing TEMP and AH values using spline
      weather_province_2019$TEMP <- na_interpolation(weather_province_2019$TEMP, option = "spline")
      weather_province_2019$AH <- na_interpolation(weather_province_2019$AH, option = "spline")

      # Save daily mean TEMP and AH for each province in 2019
      write.csv(weather_province_2019, file = paste0("data/weatherdata/weather_province_2019/", weatherdata$name[1], ".csv"), row.names = FALSE)
    }
    
  } else {
    if (!file.exists(paste0("data/weatherdata/weather_province_2019/", weatherdata$name[1], ".csv"))) {
      return(NULL)
    }
    weather_province_2019 <- read.csv(paste0("data/weatherdata/weather_province_2019/", weatherdata$name[1], ".csv"), row.names = NULL)
  }
  return(weather_province_2019)
}

# Download annual meteorological data for China from 1992-2018 and save separately
get_weather_China <- function() {
  weather_China_2002 <- get_GSOD(years = 2002, country = "China")
  save(weather_China_2002, file = "data/weatherdata/weather_China/weather_China_2002.RData")
  rm(weather_China_2002)
}
