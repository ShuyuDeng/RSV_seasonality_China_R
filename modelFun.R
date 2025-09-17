# Fit a LOESS model
callLoess <- function(dataset, variables, span) {
  loess(formula = as.formula(
          paste("AAP ~", str_c(variables, collapse = "+"), sep = "")),
        data = dataset, 
        span = span)
}

# Retrieve density value from KDE grid
calcKDE <- function(datapoint) {
  x_index <- which.min(abs(kde.grid$x - datapoint["TEMPc"]))
  y_index <- which.min(abs(kde.grid$y - datapoint["AHc"]))
  return(kde.grid$z[x_index, y_index])
}

# Predict AAP based on mean-centred temperature and absolute humidity
predictAAP <- function(TEMPc, AHc) {
  density <- calcKDE(c(TEMPc = TEMPc, AHc = AHc))
  if (density > 0) {
    predicted <- predict(Main.fit, data.frame(TEMPc = TEMPc, AHc = AHc))
    predicted[predicted<0] <- 0
    return(predicted)
  } else {
    return(NA)
  }
}

# Find the optimal threshold for predicting epidemic months
findThreshold <- function(dataset, variables = c("TEMPc", "AHc"), span = 0.3) {
  # Extract unique sites
  listOfSamples <- unique(dataset$Site)
  df <- data.frame(Site = dataset$Site,
                   Month = dataset$Month,
                   AAP = dataset$AAP,
                   Epi = dataset$Epi,
                   pred.AAP = NA,
                   matrix(data = FALSE, nrow = nrow(dataset), ncol = 19, 
                          dimnames = list(NULL, paste("Epi", 5*(1:19), sep = ""))))
  for(i in 1:length(listOfSamples)) {
    # Leave-one-site-out fitting
    fit <- callLoess(dataset[dataset$Site != listOfSamples[i],], variables, span)
    # Predict AAP for the left-out site
    predicted <- predict(fit, dataset[dataset$Site == listOfSamples[i], variables])
    predicted[predicted<0] <- 0
    pred.AAP <- predicted/sum(predicted) * 100
    df[df$Site == listOfSamples[i], 
       "pred.AAP"] <- pred.AAP
    
    # If any predicted AAP is NA, unidentified
    if(sum(is.na(pred.AAP)) >=1) {
      df[df$Site == listOfSamples[i], paste("Epi", 5*(1:19), sep = "")] <- NA
    } else {
      # Identify epidemic months based on various thresholds and cumulative predicted AAP
      for(k in 1:19){
        for (j in 1:12) {
  
          if(sum(pred.AAP[order(pred.AAP, decreasing = TRUE)[1:j]]) >= k*5) {
            df[df$Site == listOfSamples[i],][order(pred.AAP, decreasing = TRUE)[1:j], paste("Epi", k*5, sep = "")] <- TRUE
            break 
          }
        }
      }
    }
  }

  df2 <- data.frame(matrix(nrow = 19, ncol = 8))
  names(df2) <- c("threshold", "kappa", "kappa.lower", "kappa.upper", "SEN", "SPE", "PPV", "NPV")
  df2$threshold <- 5*(1:19)
  # Return Cohen's kappa, sensitivity, and specificity for each threshold
  for(k in 1:19){
    res.kappa <- psych::cohen.kappa(df[c(paste("Epi", k*5, sep = ""), "Epi")]) 
    df2$kappa[k] <- round(res.kappa$kappa, 4)
    df2$kappa.lower[k] <- round(res.kappa$confid[1,1], 4)
    df2$kappa.upper[k] <- round(res.kappa$confid[1,3], 4)

    df2$SEN[k] <- round(res.kappa$agree[2,2]/sum(res.kappa$agree[,2]) * 100, 0)
    df2$SPE[k] <- round(res.kappa$agree[1,1]/sum(res.kappa$agree[,1]) * 100, 0) 
    df2$PPV[k] <- round(res.kappa$agree[2,2] / sum(res.kappa$agree[2,]) * 100, 0) 
    df2$NPV[k] <- round(res.kappa$agree[1,1] / sum(res.kappa$agree[1,]) * 100, 0) 
  }
 return(df2)
}

# Predict epidemic months for a given site based on cumulative AAP
predictEpi <- function(site, threshold) {
  site <- site[order(site$pred.AAP, decreasing = TRUE),]
  site$cumAAP <- cumsum(site$pred.AAP)
  site <- site[order(site$Month, decreasing = FALSE),]
  site$Epi <- site$cumAAP <= threshold
  return(site)
}

# Determine RSV season onset
# It can work at different levels: city, province, or study, and supports both predicted AAP and observed AAP.
getSeasonOnset <- function(df, level) {
  
  AAP_col <- if ("pred.AAP" %in% names(df)) {
    "pred.AAP"
  } else if ("AAP" %in% names(df)) {
    "AAP"
  } else {
    stop("Neither 'pred.AAP' nor 'AAP' columns are present in the dataframe.")
  }

  if (level == "city") {
    message(df$adcode[1], df$name[1])
  } else if (level == "province") {
    message(df$Province[1])
  } else if (level == "study") {
    message(paste0(df$SID[1], ": ", df$YEAR[1]))
  } else {
    stop("Invalid level parameter. Please use 'province' or 'city'.")
  }
  
  start_row <- min(which(df$Epi == FALSE))

  if (start_row != 1) {
    df <- df %>% slice((start_row:n()) %>% c((1:(start_row-1))))
  }
  
  season_index <- 0  
  status <- 0 
  df$res <- NA
  
  for (i in 1:nrow(df)) {
    
    if (df$Epi[i] == TRUE) {
      if (status == 0) {
        season_index <- season_index + 1
        status <- 1
      }
    } else {
      status <- 0
    }
    df$res[i] <- ifelse(df$Epi[i], season_index, 0)
  }

  season_summary <- df %>%
    filter(res != 0) %>%
    group_by(res) %>%
    summarise(seasonAAP = sum(!!sym(AAP_col), na.rm = TRUE)) %>% 
    arrange(desc(seasonAAP)) %>%
    mutate(season = row_number())
  
  df <- df %>%
    left_join(season_summary, by = "res")

  add_onset <- df %>%
    filter(!is.na(season)) %>%
    group_by(season) %>%
    mutate(onset = row_number() == 1) %>%
    ungroup()

  df <- df %>%
    left_join(add_onset)  
  
  df <- df %>% arrange(Month)
  return(df)
}

# Plot RSV epidemic status maps for all months, including insets for Hong Kong and Macao
plotRSVmap <- function(df, chinamap) {
  
  month_names <- c("January", "February", "March", "April", "May", "June", 
                   "July", "August", "September", "October", "November", "December")
  
  hk_bbox <- st_bbox(chinamap %>% filter(name == "Hong Kong"))
  macao_bbox <- st_bbox(chinamap %>% filter(name == "Macao"))
  

  for (month in 1:12) {

    each_month <- df[df$Month == month,]

    map_data <- chinamap %>%
      left_join(each_month, by = c("name" = "Province", "geometry"))
    
    map_data$Status[is.na(map_data$Status)] <- "Unknown"
    
    main_map <- ggplot() +
      geom_sf(data = chinamap, color = "grey20", fill = NA, size = 5) +
      geom_sf(data = map_data, aes(fill = Status, alpha = Source), 
              color = "grey20", size = 0) +
      scale_fill_manual(values = c("Non-epidemic" = "#508AB2", 
                                   "Onset" = "#F6C14D", 
                                   "Epidemic" = "#C52A20", 
                                   "Unknown" = "white"),
                        limits = c("Onset", "Epidemic", "Non-epidemic", "Unknown"),
                        name = "Status") +
      scale_alpha_manual(values = c("Observed" = 1, 
                                    "Predicted" = 0.3),
                         na.translate = FALSE,
                         limits = c("Observed", "Predicted"),
                         name = "Source") +
      geom_hline(yintercept = 23.5, lty = 3, color = "darkgrey") +
      geom_hline(yintercept = 35, lty = 2, color = "darkgrey") +
      geom_text(aes(x = Inf, y = 23.5, label = "23.5°N"), hjust = 1, vjust = -1, size = 5, color = "grey40") +
      geom_text(aes(x = Inf, y = 35, label = "35°N"), hjust = 1, vjust = -1, size = 5, color = "grey40") +

      geom_rect(aes(xmin = hk_bbox["xmin"]-0.3, xmax = hk_bbox["xmax"]+0.1, 
                    ymin = hk_bbox["ymin"]-0.3, ymax = hk_bbox["ymax"]+0.1), 
                color = adjustcolor("red2", alpha.f = 1), fill = NA, linewidth = 0.5) +
      theme_void() +
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0)) +
      theme(legend.position = "none") +
      labs(title = month_names[month])
    
    macao_map <- ggplot() +
      geom_sf(data = map_data %>% filter(name == "Macao"), 
              aes(fill = Status, alpha = Source), 
              color = "grey20", size = 5) +
      scale_fill_manual(values = c("Non-epidemic" = "#508AB2", 
                                   "Epidemic" = "#C52A20", 
                                   "Onset" = "#F6C14D", 
                                   "Unknown" = "white"), 
                        name = "Status") +
      scale_alpha_manual(values = c("Observed" = 1, 
                                    "Predicted" = 0.3), 
                         name = "Source") +
      theme_void() +
      theme(legend.position = "none")

    hk_map <- ggplot() +
      geom_sf(data = map_data %>% filter(name == "Hong Kong"), 
              aes(fill = Status, alpha = Source), 
              color = "grey20", size = 5) +
      scale_fill_manual(values = c("Non-epidemic" = "#508AB2", 
                                   "Epidemic" = "#C52A20", 
                                   "Onset" = "#F6C14D", 
                                   "Unknown" = "white"), 
                        name = "Status") +
      scale_alpha_manual(values = c("Observed" = 1, 
                                    "Predicted" = 0.3), 
                         name = "Source") +
      theme_void() +
      theme(legend.position = "none")

    final_map <- cowplot::ggdraw() +
      draw_plot(main_map) +  
      draw_plot(macao_map, x = 0.15, y = 0.15, width = 0.1, height = 0.1) + 
      draw_plot(hk_map, x = 0.25, y = 0.15, width = 0.15, height = 0.15) +

      geom_rect(aes(xmin = 0.13, xmax = 0.42, ymin = 0.08, ymax = 0.32), 
                color = adjustcolor("red2", alpha.f = 1), fill = NA, linewidth = 0.5) +

      geom_segment(aes(x = 0.605, y = 0.372, xend = 0.13, yend = 0.32),
                   linetype = "dashed", color = "red2", alpha = 1, linewidth = 0.5) + 
      geom_segment(aes(x = 0.615, y = 0.358, xend = 0.42, yend = 0.08),
                   linetype = "dashed", color = "red2", alpha = 1, linewidth = 0.5) +

      draw_text("Macao", x = 0.19, y = 0.12, size = 14, color = "black") +
      draw_text("Hong Kong", x = 0.32, y = 0.12, size = 14, color = "black")
    

    ggsave(
      final_map,
      filename = paste("results/figures/fig4.combine_12maps/", month, ".", month_names[month], ".pdf", sep = ""),
      width = 8, height = 6
    )
  }
}

# Generate a standalone legend for RSV epidemic status maps
plotLegend <- function() {
  
  # Create dummy data for plotting legend
  dummy_data <- data.frame(
    x = 1:4,
    y = 1:4,
    Status = factor(c("Onset", "Epidemic", "Non-epidemic", "Unknown"),
                    levels = c("Onset", "Epidemic", "Non-epidemic", "Unknown")),
    Source = factor(c("Observed", "Predicted", "Observed", "Predicted"),
                    levels = c("Observed", "Predicted"))
  )
  
  # Construct legend plot
  legend_plot <- ggplot(dummy_data, aes(x = x, y = y, fill = Status, alpha = Source)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Onset" = "#F0BB41", 
                                 "Epidemic" = "#C52A20", 
                                 "Non-epidemic" = "#508AB2", 
                                 "Unknown" = "white"),
                      name = "Epidemic status: ") +
    scale_alpha_manual(values = c("Observed" = 1, "Predicted" = 0.3),
                       name = "Data source: ") +
    theme_void() + 
    theme(legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 10), 
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.spacing.x = unit(1, "cm"), 
          legend.box.just = "center",
          legend.key = element_rect(color = "black")) 
  
  # Extract the legend grob
  grobs <- ggplotGrob(legend_plot)
  legend_grob <- gtable::gtable_filter(grobs, "guide-box-bottom")
  
  pdf("results/figures/fig4.combine_12maps/legend.pdf", width = 11, height = 1.5) 
  grid.newpage()
  grid.draw(legend_grob)
  dev.off()
  
}

# Process a single StudySiteID: predict epidemic months for each study year
processStudySiteID <- function(eachStudySiteID, weather_list, studyyear) {
  
  # Filter meteorological data for the specific StudySiteID
  weather_StudySiteID <- weather_list[[eachStudySiteID]] %>%
    mutate(StudySiteID = eachStudySiteID)
  
  # Select data for each YearFrom, then combine into one dataframe
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
  
  
  # Calculate distances to study site using Vincenty ellipsoid
  weather_StudySiteID_filtered$Distance <- distVincentyEllipsoid(
    cbind(weather_StudySiteID_filtered$Longitude, weather_StudySiteID_filtered$Latitude),
    cbind(weather_StudySiteID_filtered$LONGITUDE, weather_StudySiteID_filtered$LATITUDE)
  )
  
  # Compute daily mean TEMP and AH for all weather stations within 50 km (interpolate if missing date <=5%)
  weather_StudySiteID_filtered <- weather_StudySiteID_filtered %>%
    filter(Distance <= 50000) %>%
    mutate(AH = calcAH(RH = RH, TEMP = TEMP)) %>%
    group_by(YearFrom, YEARMODA) %>%
    summarise(TEMP = mean(TEMP, na.rm = TRUE),
              AH = mean(AH, na.rm = TRUE))
  
  # Generate complete date sequence for each StudySiteID and YearFrom
  complete_dates <- studyyear %>%
    rowwise() %>%
    mutate(CompleteDate = list(seq.Date(StartDate_lag, EndDate_lag, by = "day"))) %>%
    ungroup() %>%
    select(SID, StudySiteID, Location, Province, YearFrom, CompleteDate) %>%
    unnest(CompleteDate)
  
  # Merge complete date sequence with filtered weather data and interpolate missing values
  weather_StudySiteID_complete <- complete_dates %>%
    filter(StudySiteID == eachStudySiteID) %>%
    left_join(weather_StudySiteID_filtered, by = c("YearFrom", "CompleteDate" = "YEARMODA")) %>%
    
    # Calculate missing percent, retain data with missing <=5%, and perform interpolation
    group_by(YearFrom) %>%
    mutate(MissingPercent = (sum(is.na(TEMP) | is.na(AH)) / n()) * 100) %>%
    filter(MissingPercent <= 5) %>%
    mutate(TEMP = na_interpolation(TEMP, option = "spline"),
           AH = na_interpolation(AH, option = "spline")) %>%
    select(-MissingPercent) %>%
    ungroup() %>%
    
    # Compute 14-day lagged year and month
    mutate(DatePlus14 = CompleteDate + days(14),
           Year = year(DatePlus14),
           Month = month(DatePlus14)) %>%
    select(-DatePlus14) %>% 
    
    # Compute annual average TEMP and AH for each Year, and calculate the difference between daily value and annual averages
    group_by(YearFrom) %>%
    mutate(TEMPa = mean(TEMP, na.rm = TRUE),
           AHa = mean(AH, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(TEMPc = TEMP - TEMPa, 
           AHc = AH - AHa) %>%
    
    # Compute mean-centered TEMP and AH for each year and month
    group_by(SID, StudySiteID, Location, Province, YearFrom, Year, Month) %>%
    summarise(TEMPc = mean(TEMPc, na.rm = TRUE),
              AHc = mean(AHc, na.rm = TRUE)) %>%
    ungroup() %>%
    
    # Use the predictAAP function to predict AAP for 12 months in each study year and normalize
    group_by(YearFrom) %>%
    mutate(pred.initAAP = pmap_dbl(list(TEMPc, AHc), predictAAP),
           pred.AAP = pred.initAAP / sum(pred.initAAP) * 100) %>%
    select(-pred.initAAP) %>%
    ungroup() %>%
    
    # Filter out incomplete predictions
    filter(!is.na(pred.AAP))
  
  # Determine epidemic months in each study year based on predicted AAP for 12 months and the optimal threshold (65%)
  predictEpi_results <- weather_StudySiteID_complete %>%
    group_by(SID, StudySiteID, Location, Province, YearFrom) %>%
    nest() %>%
    mutate(Epi = map(data, ~predictEpi(.x, threshold = 65))) %>%
    select(-data) %>%
    unnest(cols = c(Epi))
  
  return(predictEpi_results)
  
}

# Process all StudySiteIDs: predict epidemic months in each study year for all study sites 
processAllStudySiteIDs <- function(weather_list, studyyear, redo = FALSE) {
  
  if (redo) {
    
    # If redo, reprocesses all data; otherwise reads existing results
    StudySiteIDs <- names(weather_list)
    
    results <- lapply(StudySiteIDs, function(eachStudySiteID) {
      processStudySiteID(eachStudySiteID, weather_list, studyyear)
    })
    
    final_result <- bind_rows(results)
    write.csv(final_result, "results/tables/predictEpi_study_results.csv", row.names = FALSE)
    
    return(final_result)
  } else {
    
    # Read previously saved results
    return(read.csv("results/tables/predictEpi_study_results.csv"))
  }
}


