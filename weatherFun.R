library(GSODR)
library(geosphere)
library(imputeTS)

# 获取距离每个研究地点200公里内气象站的所有研究期限内的所有气象数据
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
calcAH <- function(RH, TEMP) {
  AH <- 6.112*exp((17.67*TEMP)/(TEMP + 243.5))*RH*2.1674/(273.15 + TEMP)
  return(AH)
}

# 对每个研究地点的气象数据进行清理，如果有多个年份，计算每月每日的均值
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

get_weather_China_2019 <- function() {
  # 中国2019气象数据：缺澳门
  weather_China_2019 <- get_GSOD(years = 2019, country = "China")
  unique(weather_China_2019$Province)
  save(weather_China_2019, file = "data/weatherdata/weather_China/weather_China_2019.RData")
  # 补充澳门2019气象数据
  weather_Macao_2019 <- get_GSOD(years = 2019, country = "MO")
  save(weather_Macao_2019, file = "data/weatherdata/weather_China/weather_Macao_2019.RData")
  rm(weather_China_2019, weather_Macao_2019)
}

get_weather_province_2019 <- function(weatherdata, redo = FALSE) {
  message(weatherdata$Province[1], weatherdata$name[1])

  if (redo) {
    # 计算绝对湿度，按省份分组，计算2019年每日气象均值
    weather_province_2019 <- weatherdata %>%
      mutate(AH = calcAH(RH = RH, TEMP = TEMP)) %>%
      group_by(YEARMODA) %>%
      summarise(Province = Province[1],
                DaySeq = YDAY[1],
                TEMP = mean(TEMP, na.rm = TRUE), 
                AH = mean(AH, na.rm = TRUE)) %>%
      mutate(DaySeq = order(YEARMODA))
    missing.percent <- (365 - n_distinct(weather_province_2019$DaySeq))/365 * 100
    # 暂定：如果缺失比例超过5%，不返回该省份2019气象数据
    if (missing.percent > 5) {
      return(NULL)
    } else {
      # 否则创建一个完整的日期数据框，添加省份信息
      completeDates <- seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "day")
      df <- data.frame(Province = weatherdata$Province[1],
                       YEARMODA = completeDates, 
                       DaySeq = as.numeric(format(completeDates, "%j")))
      weather_province_2019 <- df %>%
        left_join(weather_province_2019, by = c("Province", "YEARMODA", "DaySeq"))
      # 对温湿度进行spline插补
      weather_province_2019$TEMP <- na_interpolation(weather_province_2019$TEMP, option = "spline")
      weather_province_2019$AH <- na_interpolation(weather_province_2019$AH, option = "spline")
      # 保存各省份2019年每日平均温湿度数据
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


# 下载中国1992-2018年的每年气象数据并分别保存
get_weather_China <- function() {
  
  weather_China_2002 <- get_GSOD(years = 2002, country = "China")
  save(weather_China_2002, file = "data/weatherdata/weather_China/weather_China_2002.RData")
  rm(weather_China_2002)
  
}

# 处理单个StudySiteID
processStudySiteID <- function(eachStudySiteID, weather_list, studyyear) {
  # 筛选并标记每个研究完整时间段的气象数据
  weather_StudySiteID <- weather_list[[eachStudySiteID]] %>%
    mutate(StudySiteID = eachStudySiteID)
  
  # 对每个YearFrom进行筛选并标记，然后将结果合并为一个数据框
  weather_StudySiteID_filtered <- studyyear %>%
    filter(StudySiteID == eachStudySiteID) %>%
    rowwise() %>%
    do({
      data <- filter(weather_StudySiteID, YEARMODA >= .$StartDate_lag & YEARMODA <= .$EndDate_lag)
      mutate(data, YearFrom = .$YearFrom,
             Latitude = .$Latitude, 
             Longitude = .$Longitude)
    }) %>%
    ungroup() %>%
    bind_rows()
  
  
  # 向量化计算距离
  weather_StudySiteID_filtered$Distance <- distVincentyEllipsoid(
    cbind(weather_StudySiteID_filtered$Longitude, weather_StudySiteID_filtered$Latitude),
    cbind(weather_StudySiteID_filtered$LONGITUDE, weather_StudySiteID_filtered$LATITUDE)
  )
  
  # 计算距离每个研究50公里内所有气象站，对应时间段的日均TEMP+AH（如缺失天数少于5%，进行插补）
  weather_StudySiteID_filtered <- weather_StudySiteID_filtered %>%
    filter(Distance <= 50000) %>%
    mutate(AH = calcAH(RH = RH, TEMP = TEMP)) %>%
    group_by(YearFrom, YEARMODA) %>%
    summarise(TEMP = mean(TEMP, na.rm = TRUE),
              AH = mean(AH, na.rm = TRUE))
  
  # 生成每个StudySiteID和YearFrom的完整日期序列
  complete_dates <- studyyear %>%
    rowwise() %>%
    mutate(CompleteDate = list(seq.Date(StartDate_lag, EndDate_lag, by = "day"))) %>%
    ungroup() %>%
    select(SID, StudySiteID, Location, Province, YearFrom, CompleteDate) %>%
    unnest(CompleteDate)
  
  # 合并日期序列和筛选后的数据
  weather_StudySiteID_complete <- complete_dates %>%
    filter(StudySiteID == eachStudySiteID) %>%
    left_join(weather_StudySiteID_filtered, by = c("YearFrom", "CompleteDate" = "YEARMODA")) %>%
    
    # 计算缺失比例，保留缺失比例<=5%的数据，进行TEMP和AH的插补
    group_by(YearFrom) %>%
    mutate(MissingPercent = (sum(is.na(TEMP) | is.na(AH)) / n()) * 100) %>%
    filter(MissingPercent <= 5) %>%
    mutate(TEMP = na_interpolation(TEMP, option = "spline"),
           AH = na_interpolation(AH, option = "spline")) %>%
    select(-MissingPercent) %>%    # 移除缺失百分比辅助列
    ungroup() %>%
    
    # 计算滞后14天日期对应的月份，新增Year和Month列
    mutate(DatePlus14 = CompleteDate + days(14),
           Year = year(DatePlus14),
           Month = month(DatePlus14)) %>%
    select(-DatePlus14) %>%  # 移除临时列DatePlus14
    
    # 计算每个YearFrom的年均TEMP和AH, 并计算每日中心差值
    group_by(YearFrom) %>%
    mutate(TEMPa = mean(TEMP, na.rm = TRUE),
           AHa = mean(AH, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(TEMPc = TEMP - TEMPa, 
           AHc = AH - AHa) %>%
    
    # 计算每个研究年限中，每月的中心差值均值
    group_by(SID, StudySiteID, Location, Province, YearFrom, Year, Month) %>%
    summarise(TEMPc = mean(TEMPc, na.rm = TRUE),
              AHc = mean(AHc, na.rm = TRUE)) %>%
    ungroup() %>%
    
    # 使用predictAAP函数预测每个研究年限中，12个月的AAP，并直接进行归一化处理
    group_by(YearFrom) %>%
    mutate(pred.initAAP = pmap_dbl(list(TEMPc, AHc), predictAAP),
           pred.AAP = pred.initAAP / sum(pred.initAAP) * 100) %>%
    select(-pred.initAAP) %>%
    ungroup() %>%
    
    # 筛选出可预测12个月AAP的研究年限
    filter(!is.na(pred.AAP))
  
  # 基于12个月AAP和预测模型的最佳阈值65%，判断各研究每个研究年限的epidemic months
  predictEpi_results <- weather_StudySiteID_complete %>%
    group_by(SID, StudySiteID, Location, Province, YearFrom) %>%
    nest() %>%
    mutate(Epi = map(data, ~predictEpi(.x, threshold = 65))) %>%
    select(-data) %>%
    unnest(cols = c(Epi))
  
  return(predictEpi_results)
  
}

processAllStudySiteIDs <- function(weather_list, studyyear, redo = FALSE) {
  
  if (redo) {
    
    # 重新处理数据并保存结果
    StudySiteIDs <- names(weather_list)
    
    results <- lapply(StudySiteIDs, function(eachStudySiteID) {
      processStudySiteID(eachStudySiteID, weather_list, studyyear)
    })
    
    final_result <- bind_rows(results)
    write.csv(final_result, "results/tables/predictEpi_study_results.csv", row.names = FALSE)
    
    return(final_result)
  } else {
    
    # 读取已有CSV文件
    return(read.csv("results/tables/predictEpi_study_results.csv"))
    }
}
