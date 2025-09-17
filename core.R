rm(list=ls())
library(stringr)
library(ggplot2)
library(sf)
library(jsonlite)
library(MASS)
library(RcppRoll) 
library(cowplot)
library(ggpattern)
source("dataProcessing.R")
source("modelFun.R")

# 0. Default----
density_threshold <- 0.0035
bestLag <- 14

# 1. Description----
# 1.1 Study quality----
# No. of low-quality studies：21
metadata %>% 
  filter(QA1 == "C" | QA2 == "C" |QA3 == "C") %>% 
  summarise(n_SID = n_distinct(SID))

# 1.2 Sample size----
# Only sites reporting >=100 RSV cases included
table_summary_case2 <- RSVdata3 %>%
  group_by(Location, Province, Latitude, Longitude) %>%
  summarise(TotalRSV = sum(SumRSV))
# No. of RSV cases：193843
sum(table_summary_case2$TotalRSV)
# No. of study sites：44
nrow(table_summary_case2)
# No. of provinces：21
length(unique(table_summary_case2$Province))
# No. of study-years: 427 (for each site: 1-42)
table_summary_year2 <- RSVdata4 %>%
  filter(Month == "1") %>%
  group_by(Location, Province, Latitude, Longitude) %>%
  summarise(TotalYear = sum(SumYear))
sum(table_summary_year2$TotalYear)
range(table_summary_year2$TotalYear)

# 1.3 Map of study sites----
chinamap_ADM1 <- sf::st_read("data/chinamap/ADM1/中华人民共和国.shp")
chinamap_ADM1$name <- sub("Nei Mongol", "Inner Mongolia", chinamap_ADM1$name)

# Convert study sites to sf
sites_sf <- sf::st_as_sf(RSVdata4 %>% distinct(Latitude, Longitude, Province, SumYear),
                     coords = c("Longitude", "Latitude"), 
                     crs = st_crs(chinamap_ADM1))

# Mark provinces with study sites
chinamap_ADM1$withsite <- as.logical(chinamap_ADM1$name %in% sites_sf$Province)

# Get bbox for Hong Kong & Macao
hk_bbox <- sf::st_bbox(chinamap_ADM1 %>% filter(name == "Hong Kong"))
macao_bbox <- sf::st_bbox(chinamap_ADM1 %>% filter(name == "Macao"))

# Main map
main_map <- ggplot() +
  geom_sf(data = chinamap_ADM1, color = "darkgrey", aes(fill = withsite)) +
  scale_fill_manual(values = c("FALSE" = "gainsboro", "TRUE" = "aliceblue"),
                    name = "Data availability",
                    labels = c("Not available", "Available")) +
  geom_sf(data = sites_sf, aes(size = SumYear), color = "red3", alpha = 0.5) +
  scale_size_continuous(range = c(1, 3), name = "Years of data") +
  geom_sf_text(data = chinamap_ADM1 %>% drop_na(name), aes(label = ifelse(withsite, name, "")), size = 2) +
  geom_hline(yintercept = 23.5, lty = 3, color = "darkgrey") +
  geom_hline(yintercept = 35, lty = 2, color = "darkgrey") +
  geom_text(aes(x = Inf, y = 23.5, label = "23.5°N"), hjust = 1, vjust = -1, size = 3, color = "grey40") +
  geom_text(aes(x = Inf, y = 35, label = "35°N"), hjust = 1, vjust = -1, size = 3, color = "grey40") +
  
  # Highlight Hong Kong & Macao with red boxes
  geom_rect(aes(xmin = hk_bbox["xmin"]-0.3, xmax = hk_bbox["xmax"]+0.1, 
                ymin = hk_bbox["ymin"]-0.3, ymax = hk_bbox["ymax"]+0.1), 
            color = adjustcolor("red2", alpha.f = 1), fill = NA, linewidth = 0.5) +
  theme_void() + 
  theme(legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 8),
        legend.position = c(0.9, 0.2),
        legend.box = "horizontal")

# Insets: Hong Kong & Macao
hk_map <- ggplot() +
  geom_sf(data = chinamap_ADM1 %>% filter(name == "Hong Kong"), color = "darkgrey", aes(fill = withsite)) +
  scale_fill_manual(values = c("FALSE" = "gainsboro", "TRUE" = "aliceblue"),
                    name = "Data availability") +
  geom_sf(data = sites_sf %>% filter(Province == "Hong Kong"), aes(size = SumYear), color = "red3", alpha = 0.5) +
  theme_void() +
  theme(legend.position = "none")
macao_map <- ggplot() +
  geom_sf(data = chinamap_ADM1 %>% filter(name == "Macao"), color = "darkgrey", aes(fill = withsite)) +
  scale_fill_manual(values = c("FALSE" = "gainsboro", "TRUE" = "aliceblue"),
                    name = "Data availability") +
  geom_sf(data = sites_sf %>% filter(Province == "Macao"), aes(size = SumYear), color = "red3", alpha = 0.5) +
  theme_void() +
  theme(legend.position = "none")

# Combine main map and insets
final_map <- cowplot::ggdraw() +
  draw_plot(main_map) +
  # Red frame for inset area
  geom_rect(aes(xmin = 0.15, xmax = 0.42, ymin = 0.08, ymax = 0.32), 
            color = "red2", fill = "transparent", linewidth = 0.5) +
  # Insets
  draw_plot(macao_map, x = 0.15, y = 0.15, width = 0.1, height = 0.1) +  
  draw_plot(hk_map, x = 0.25, y = 0.15, width = 0.15, height = 0.15) + 
  # Connection lines
  geom_segment(aes(x = 0.61, y = 0.394, xend = 0.15, yend = 0.32),
               linetype = "dashed", color = "red2", alpha = 1, linewidth = 0.5) +  
  geom_segment(aes(x = 0.622, y = 0.38, xend = 0.42, yend = 0.08),
               linetype = "dashed", color = "red2", alpha = 1, linewidth = 0.5) +  
  # Text labels
  draw_text("Macao", x = 0.2, y = 0.12, size = 10, color = "black") +  
  draw_text("Hong Kong", x = 0.32, y = 0.12, size = 10, color = "black")  
ggsave("results/figures/fig1.study_sites.pdf", final_map, width = 8, height = 6,  bg = "white")

# 2. Modelling----
# 2.1 Model selection----
# Re-run "compareModels" for model selection
# Best model: TEMPc + AHc, span = 0.30, lag = 14 days

# 2.2 Model fit----
# Scatter plot of TEMPc and AHc: predictors are correlated and occupy only a limited region
plot(x = dataFinal$Lag_14$TEMPc, y = dataFinal$Lag_14$AHc,
     xlab = "TEMPc (°C)", ylab = "AHc (g/m³)")

# use KDE to define a data-supported prediction domain
# Rationale: LOESS can behave poorly when asked to predict in regions with little/no data support. 
# We therefore estimate the joint (TEMPc, AHc) density and restrict predictions to areas with sufficient support.
# 2D Gaussian KDE on a 100×100 grid
kde.grid <- with(dataFinal$Lag_14, 
                 MASS::kde2d(TEMPc, AHc, n = 100))
contour(kde.grid) # Diagnostic plot of the KDE surface

# Model fit
Main.fit <- callLoess(dataFinal$Lag_14, c("TEMPc", "AHc"), span = 0.3)

# Build a fine prediction grid and attach a KDE value to each (TEMPc, AHc) point
Main.fit.grid <- expand.grid(
  TEMPc = seq(-21, 17, 0.1), # Set ranges to cover the observed domain (± small margins)
  AHc = seq(-8, 12, 0.1)
) %>%
  rowwise() %>%
  mutate(density = calcKDE(c(TEMPc = TEMPc, AHc = AHc))) %>%  # KDE support at this point
  filter(density > density_threshold) %>%                     # Keep only well-supported points
  ungroup() %>% 
  mutate(AAP = predict(Main.fit, newdata = dplyr::select(., TEMPc, AHc)))   # predict AAP
summary(Main.fit.grid)
Main.fit.grid$AAP[Main.fit.grid$AAP < 0] <- 0

# Figure 2
ggsave(
  ggplot(data = Main.fit.grid, aes(TEMPc, AHc)) +
    geom_point(aes(color = AAP)) +
    scale_color_gradient2(low = "dodgerblue3", mid = "white", high = "red3",
                          limits = c(0, 30), midpoint = 8,
                          name = "AAP (%)") +
    scale_x_continuous(breaks = c(-24, -16, -8, 0, 8, 16, 24), expand = c(0,0)) +
    scale_y_continuous(breaks = c( -12, -8, -4, 0, 4, 8, 12), expand = c(0,0)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = "Mean-centred temperature (°C)", 
         y = expression("Mean-centred absolute humidity (g/m"^3*")")) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12),  # X-axis tick labels
          axis.text.y = element_text(size = 12),  
          axis.title.x = element_text(size = 12), # X-axis title
          axis.title.y = element_text(size = 12))
  ,
  filename = "results/figures/fig2.model_fit.pdf", width = 6, height = 5
)

# 2.3 Cut-off selection----
# 1. Calculate kappa, sensitivity, specificity, PPV, and NPV for 19 candidate thresholds
thresholds <- findThreshold(dataFinal$Lag_14)

# 2. Plot ROC curve for the 19 thresholds and mark the optimal cut-off point
# Add points (0,0) and (100,100) 
thresholds <- rbind(data.frame(threshold = 0, 
                               kappa = NA, kappa.lower = NA, kappa.upper = NA, 
                               SEN = 0, SPE = 100, PPV = NA, NPV = NA), 
                    thresholds)
thresholds <- rbind(thresholds, 
                    data.frame(threshold=100, 
                               kappa = NA, kappa.lower = NA, kappa.upper = NA, 
                               SEN = 100, SPE = 0, PPV = NA, NPV = NA))

# Convert sensitivity and specificity into ROC coordinates
thresholds$TPR <- thresholds$SEN         # True Positive Rate for the y-axis
thresholds$FPR <- 100 - thresholds$SPE   # False Positive Rate for the x-axis

# Sort thresholds in ascending order of FPR
thresholds <- thresholds[order(thresholds$FPR),]

# Calculate AUC using the trapezoidal rule
# Note: FPR and TPR are in percentages, so divide by 100*100
auc <- sum((thresholds$FPR[-1] - thresholds$FPR[-length(thresholds$FPR)]) * 
             (thresholds$TPR[-1] + thresholds$TPR[-length(thresholds$TPR)]) / 2) / (100 * 100)

# Figure S1: ROC curve
ggsave(
  ggplot(thresholds, aes(x = FPR, y = TPR)) +
    geom_line(color = "dodgerblue2") +
    geom_point(color = "dodgerblue4", alpha = 0.6, size=2) +
    
    # Mark the optimal cut-off point (maximum kappa) in red
    geom_point(data = thresholds[which.max(thresholds$kappa),],   
               aes(x = FPR, y = TPR), color="red", alpha = 0.6, size=2) +
    
    # Add diagonal reference line
    geom_segment(aes(x = 0, y = 0, xend = 100, yend = 100), linetype = "dashed", color = "grey") +
    labs(x = "1 - Specificity (False Positive Rate, %)", 
         y = "Sensitivity (True Positive Rate, %)") +
    
    # Add annotation for the optimal threshold value (65%)
    annotate("text", x=thresholds$FPR[which.max(thresholds$kappa)], y=thresholds$TPR[which.max(thresholds$kappa)], 
             label=paste(thresholds$threshold[which.max(thresholds$kappa)], "%", sep = ""), 
             hjust=1, vjust=-1.4, color = "red", size = 5) +
    
    # Add annotation for the AUC value
    annotate("text", x = 75, y = max(thresholds$TPR)/3, 
             label = paste("AUC =", round(auc, 4)), size = 7, color = "dodgerblue4") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12))
  ,
  filename = "results/figures/figs1.ROC_for_19thresholds.pdf", width = 6, height = 6
)

# 3. Prediction at province-level----
# 3.1 Predict AAP by province-month----
# Check the coverage of Chinese meteorological data in 2019: 34 provinces
unique(weather_China_2019$Province)

# Aggregate daily meteorological data by province
#  - Calculate daily mean TEMP and AH across all weather stations within each province
#  - If missing days < 5%, perform imputation
weather_province_2019 <- do.call(rbind, by(weather_China_2019, weather_China_2019$Province, get_weather_province_2019,
                                           redo = FALSE))
# Check for unexpected missing province values: none
weather_province_2019[which(is.na(weather_province_2019$Province)),]

# Compute mean-centred TEMP (TEMPc) and mean-centred AH (AHc), and predict AAP for each province-month
predict_province <- weather_province_2019 %>%
  # Step 1. Calculate province-level annual means of TEMP and AH
  group_by(Province) %>%
  summarise(TEMPa = mean(TEMP, na.rm = TRUE),
            AHa = mean(AH, na.rm = TRUE)) %>%
  
  # Step 2. Merge with original daily data and compute daily TEMPc and AHc
  right_join(weather_province_2019, by = "Province") %>%
  mutate(TEMPc = TEMP - TEMPa, 
         AHc = AH - AHa,
         Month = cut(ifelse(DaySeq - bestLag < 1, DaySeq - bestLag + 365, DaySeq - bestLag),
                     breaks = c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365),
                     labels = c(1:12), right = TRUE)) %>%
  
  # Step 3. Calculate lagged monthly averages of TEMPc and AHc
  group_by(Province, Month) %>%
  summarise(TEMPc = mean(TEMPc, na.rm = TRUE), 
            AHc = mean(AHc, na.rm = TRUE)) %>%

  # Step 4. Predict AAP for each province-month
  mutate(pred.initAAP = pmap_dbl(list(TEMPc, AHc), predictAAP),
         # Normalise to ensure AAP sums to 100% across 12 months
         pred.AAP = pred.initAAP / sum(pred.initAAP) * 100) %>%
  dplyr::select(-pred.initAAP) %>%
  ungroup() %>%
  
  # Step 5. Retain provinces with complete 12-month predictions
  filter(!is.na(pred.AAP))

# Check predicted provinces: 31 provinces in total 
# 21 provinces with observed data + 10 provinces with predictions only
unique(predict_province$Province)

# 3.2 Identify predicted epidemic months and onset by province----
# Epidemic months
#  - Based on predicted monthly AAP values
#  - Using the optimal threshold (65%) from the prediction model
predict_province2 <- do.call(rbind, 
                             by(predict_province, predict_province$Province, predictEpi, threshold = 65))
rownames(predict_province2) <- NULL

# Season onset
res.province <- do.call(rbind, 
                        by(predict_province2, predict_province2$Province, getSeasonOnset, level = "province"))

# Classify epidemic status
res.province$Status <- ifelse(is.na(res.province$season), "Non-epidemic", 
                              ifelse(res.province$onset, "Onset", "Epidemic"))
# Merge results with province-level shapefile (ADM1) for mapping
res.province <- left_join(res.province, chinamap_ADM1, by = c("Province" = "name"))
summary(res.province)

# 3.3 Comparison of observed and predicted epidemic timing across 21 provinces arranged by latitude----
# Clean observed data
# Extract province centroids from ADM2 and add to ADM1.1
chinamap_ADM2 <- list.files(path = "data/chinamap/ADM2", pattern = ".json", full.names = TRUE) %>%
  map_dfr(~ st_read(.x) %>% 
            mutate(adcode = as.character(adcode)))
filtered_data <- chinamap_ADM2 %>%
  filter(level == "province") %>%
  head(-1) %>%   # drop last row (duplicated)
  dplyr::select(center) %>%
  unique()
chinamap_ADM1.1 <- chinamap_ADM1 %>%
  mutate(Longitude = c(sapply(filtered_data$center, function(x) x[1]), NA),
         Latitude = c(sapply(filtered_data$center, function(x) x[2]), NA))

# Process published province-level observed RSV data by Guo et al.
RSVprovincedata2 <- RSVprovincedata_pub %>%
  # Keep provinces with at least 100 RSV detections (remove Ningxia and Tianjin)
  filter(RSVn >= n.sample) %>%
  dplyr::select(Province, starts_with("m")) %>%
  # Reshape to long format (Month, Status)
  pivot_longer(cols = starts_with("m"), 
               names_to = "Month", 
               values_to = "Status") %>%
  mutate(Month = as.integer(gsub("m", "", Month))) %>%
  # Add coordinates and sort by latitude
  left_join(chinamap_ADM1.1[c("name", "Longitude", "Latitude")], 
            by = c("Province" = "name")) %>%
  arrange(Latitude)
# Keep only epidemic and onset months for comparison
observed_data <- RSVprovincedata2 %>% 
  dplyr::select(Province, Month, Status, Longitude, Latitude)%>%
  filter(Status %in% c("Epidemic", "Onset"))
province_compare <- unique(observed_data$Province)

# Clean predicted data
predicted_data <- res.province %>%
  # Add coordinates and sort by latitude
  left_join(chinamap_ADM1.1[c("name", "Longitude", "Latitude", "geometry")],
            by = c("Province" = "name", "geometry")) %>%
  arrange(Latitude) %>%
  # Keep only predicted epidemic and onset months
  filter(Province %in% province_compare) %>% 
  dplyr::select(Province, Month, Status, Longitude, Latitude) %>%
  filter(Status %in% c("Epidemic", "Onset"))

# Combine observed and predicted data
# Distinguish data source: observed vs predicted
observed_data$Source = 'Observed'
predicted_data$Source = 'Predicted'
# Map provinces to numeric positions for y-axis alignment
province_mapping <- data.frame(Province = province_compare, 
                               Num = 1:length(province_compare))
# Assign Y-axis positions (Observed above, Predicted below)
observed_data <- observed_data %>% 
  left_join(province_mapping, by = "Province") %>%
  mutate(Province_Num = Num + 0.1)
predicted_data <- predicted_data %>% 
  left_join(province_mapping, by = "Province") %>%
  mutate(Province_Num = Num - 0.1)
combined_data <- rbind(observed_data, predicted_data)

# Count epidemic months per source-province
combined_data <- combined_data %>%
  group_by(Source, Province) %>%
  mutate(EpiMonthsCount = n()) %>%
  ungroup()

# Mark provinces with clear seasonality (with ≤5 epidemic months in observed data) with an asterisk
combined_data <- combined_data %>%
  mutate(ProvinceStar = ifelse(Source == "Observed" & EpiMonthsCount <= 5,
                               paste0(Province, "*"), Province))

# Create unique labels for y-axis (ProvinceStar + Latitude, keep order)
combined_data$ProvinceLat <- paste0(combined_data$ProvinceStar, "\n(", round(combined_data$Latitude, 1), ")")
provinceLat_compare <- combined_data %>%
  filter(Source == "Observed") %>%
  pull(ProvinceLat) %>%
  unique()

# Create legend (Source + Status)
combined_data$SourceStatus <- paste(combined_data$Source, combined_data$Status, sep = " ")
combined_data$SourceStatus <- factor(combined_data$SourceStatus, 
                                   levels = c("Observed Onset", "Observed Epidemic", "Predicted Onset", "Predicted Epidemic"))
# Define color scheme
color_mapping <- c(
  "Observed Onset" = "#92B5CA",
  "Observed Epidemic" = "#E69191",
  "Predicted Onset" = "#CCE4EF",
  "Predicted Epidemic" = "#F5DFDB"
)

# Figure 3
ggsave(
  ggplot(combined_data, aes(x = as.numeric(Month)-0.5, xend = as.numeric(Month)+0.5, 
                            y = Province_Num, yend = Province_Num, color = SourceStatus)) +
    geom_segment(linewidth = 2.5) +
    scale_color_manual(values = color_mapping, name = "Status") +
    scale_x_continuous(breaks = c(1:12), name = "Month") +
    scale_y_continuous(breaks = province_mapping$Num, labels = provinceLat_compare, name = "Province (Latitude)") +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(), 
          panel.grid.minor = element_blank())
  ,
  filename = "results/figures/fig3.compare_onset_epi.pdf", width = 6, height = 8
)

# 3.4 Observed plus predicted RSV epidemic timing across 31 provinces----
combined_map <- rbind(
  # 21 provinces with observed data
  RSVprovincedata2 %>%
    dplyr::select(Province, Month, Status, geometry) %>%    
    mutate(Source = "Observed"),
  # 10 provinces with predicted data
  res.province %>%
    anti_join(RSVprovincedata2, by = "Province") %>%
    dplyr::select(Province, Month, Status, geometry) %>%    
    mutate(Source = "Predicted")
  )

# Figure 4
plotRSVmap(combined_map, chinamap_ADM1)

# Filter provinces had RSV season onset between October and December: 26/31, 84%
combined_map %>% 
  filter(Status == "Onset", Month %in% c("10", "11", "12")) %>%
  pull(Province) %>%
  unique()

# 4. Exploratory analysis: comparison of predicted and observed year-on-year variations----
# 4.1 Observed year-on-year variations----
# Select study sites with RSV activity data for ≥3 consecutive years: RSVdata_2
# Identify epidemic months for each site-year
RSVdata_3 <- do.call(rbind, 
                     by(RSVdata_2, RSVdata_2[, c("StudySiteID", "YearFrom")], getAAPEpi, RSVn = "RSV", threshold = 75))

# Create a date column (Year-Month-01) for plotting x-axis, sorted by study site and date
RSVdata_3 <- RSVdata_3 %>%
  mutate(YearMonth = as.Date(paste(Year, Month, "01", sep = "-"))) %>%
  arrange(StudySiteID, YearMonth) 

# 4.2 Predicted year-on-year variations----
# Prepare 14-day lagged weather data
RSVdata_4 <- RSVdata_3 %>%
  mutate(StartDate_lag = as.Date(paste(YearFrom, MonthFrom, "01", sep = "-")) - days(14),
         EndDate_lag = (as.Date(paste(YearTo, MonthTo, "01", sep = "-")) %>% 
                          ceiling_date(unit = "month") - days(1)) - days(14)
  )

# Determine required years of weather data: 1992–2019
sort(unique(RSVdata_4$StartDate_lag))
sort(unique(RSVdata_4$EndDate_lag))

# Weather data extraction: within 50 km of study site
# Load pre-downloaded China weather data for 1992–2019
file_list <- list.files("data/weatherdata/weather_China", full.names = TRUE)
for (file in file_list) {
  load(file)
}
weather_China <- bind_rows(
  mget(ls(pattern = "^weather_China_\\d{4}$"))
)

# Prepare a list to store site-specific weather data for each study site
weather_list <- list()
for(each in unique(RSVdata_4$StudySiteID)) {
  time_range <- RSVdata_4 %>%
    filter(StudySiteID == each) %>%
    summarise(Start = min(StartDate_lag), End = max(EndDate_lag))
  weather_subset <- weather_China %>%
    filter(YEARMODA >= time_range$Start & YEARMODA <= time_range$End)
  weather_list[[each]] <- weather_subset
}

# Predict epidemic months for each study site and year
studyyear <- RSVdata_4 %>%
  group_by(SID, StudySiteID, Latitude, Longitude, Location, Province, YearFrom) %>%
  summarise(
    StartDate_lag = unique(StartDate_lag),
    EndDate_lag = unique(EndDate_lag),
  )
predictEpi_study <- processAllStudySiteIDs(weather_list, studyyear, redo = FALSE)

# 4.3 Comparison----
# Identify study sites with ≥3 consecutive years of predicted results
consecutive_yearfrom <- predictEpi_study %>%
  group_by(StudySiteID) %>%
  summarise(UniqueYearFrom = unique(YearFrom)) %>%
  arrange(StudySiteID, UniqueYearFrom) %>%
  mutate(ConsecutiveDiff = c(1, diff(UniqueYearFrom)),      # year-to-year gap
         Group = cumsum(ConsecutiveDiff != 1),              # new group if not consecutive
         ConsecutiveCount = ave(UniqueYearFrom, StudySiteID, Group, 
                                FUN = length)) %>%
  
  filter(ConsecutiveCount >= 3) %>%    # keep ≥3 consecutive years
  ungroup() %>%
  dplyr::select(-ConsecutiveDiff, -Group, -ConsecutiveCount,
                YearFrom = UniqueYearFrom) 

# 17 studies met the criterion
length(unique(consecutive_yearfrom$StudySiteID))

# Filter predicted epidemic results to keep only the consecutive years
predictEpi_study_consecutive <- predictEpi_study %>%
  semi_join(consecutive_yearfrom, by = c("StudySiteID", "YearFrom")) %>%
  mutate(YearMonth = as.Date(paste(Year, Month, "01", sep = "-"))) %>%
  arrange(StudySiteID, YearFrom, YearMonth)  %>% 
  dplyr::select(SID, StudySiteID, Location, Province, YearFrom, YearMonth, Epi) %>%
  mutate(Source = "Predicted")

# Filter observed epidemic results in the same way
RSVdata_3_consecutive <- RSVdata_3 %>%
  semi_join(consecutive_yearfrom, by = c("StudySiteID", "YearFrom")) %>%
  arrange(StudySiteID, YearFrom, YearMonth) %>% 
  dplyr::select(SID, StudySiteID, Location, Province, YearFrom, YearMonth, Epi) %>%
  mutate(Source = "Observed")

# Combine observed and predicted results
combined_consecutive <- rbind(RSVdata_3_consecutive, predictEpi_study_consecutive) %>%
  mutate(Status = ifelse(Epi == TRUE, "Epidemic", "Non-epidemic")) %>%
  mutate(SourceStatus = paste(Source, Status, sep = " ")) %>%
  left_join(metadata[c("SID", "Author", "PubYear")])

# Prepare for plotting
combined_consecutive <- combined_consecutive %>%
  mutate(
    SourceLevel = ifelse(Source == "Predicted", 0.3, 0.7),  # assign Y-axis positions (Observed above, Predicted below)
    MonthStart = YearMonth,
    MonthEnd = MonthStart %m+% months(1) ,  # ensure continuous monthly segments
    FacetLabel = ifelse(
      is.na(Location),
      paste(Province, " (", Author, " ", PubYear, ")", sep = ""),
      paste(Location, ", ", Province, " (", Author, " ", PubYear, ")", sep = "")
    )
  )
unique(combined_consecutive[, c("SID", "Author", "FacetLabel")])

#  Identify study sites with clear seasonality in observed data
seasonality_consecutive <- combined_consecutive %>%
  filter(Source == "Observed") %>%
  group_by(StudySiteID, YearFrom) %>%
  summarize(EpiCount = sum(Epi == TRUE, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(StudySiteID) %>%
  mutate(ClearSeasonality = all(EpiCount <= 5)) %>%  # clear seasonality
  ungroup() %>%
  dplyr::select(-EpiCount) %>% 
  left_join(combined_consecutive, by = c("StudySiteID", "YearFrom"))

# Define color scheme
color_mapping <- c(
  "Observed Non-epidemic" = "#92B5CA",
  "Observed Epidemic" = "#E69191",
  "Predicted Non-epidemic" = "#CCE4EF",
  "Predicted Epidemic" = "#E6CECF"
)

# Figure S2. (A) study sites with clear seasonality
fig_a <- ggplot(subset(seasonality_consecutive, ClearSeasonality == TRUE)) +
  geom_segment(aes(x = MonthStart, xend = MonthEnd,
                   y = SourceLevel, yend = SourceLevel, color = SourceStatus), size = 4) +
  facet_wrap(~FacetLabel, scales = "free_x", ncol = 1) +
  labs(x = "Year-month", y = "") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y-%m", minor_breaks = "1 month") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0.3, 0.7), 
                     labels = c("Predicted", "Observed")) +
  scale_color_manual(values = color_mapping, name = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 13),
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        panel.grid.major.x = element_line(color="grey", size=0.5),
        panel.grid.minor.x = element_line(color="lightgrey", size=0.25),
        legend.position = "bottom",
        plot.margin = margin(t = 20, r = 20, unit = "pt")
  )

# Figure S2. (B) study sites with no clear seasonality
fig_b <- ggplot(subset(seasonality_consecutive, ClearSeasonality == FALSE)) +
  geom_segment(aes(x = MonthStart, xend = MonthEnd,
                   y = SourceLevel, yend = SourceLevel, color = SourceStatus), size = 4) +
  facet_wrap(~FacetLabel, scales = "free_x", ncol = 1) +
  labs(x = "Year-month", y = "") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y-%m", minor_breaks = "1 month") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0.3, 0.7), 
                     labels = c("Predicted", "Observed")) +
  scale_color_manual(values = color_mapping, name = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 13),
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        panel.grid.major.x = element_line(color="grey", size=0.5),
        panel.grid.minor.x = element_line(color="lightgrey", size=0.25),
        legend.position = "bottom",
        plot.margin = margin(t = 20, r = 20, unit = "pt")
  )

# Remove raw legends and extract legend separately
fig_a_nolegend <- fig_a + theme(legend.position = "none")
fig_b_nolegend <- fig_b + theme(legend.position = "none")
grob <- ggplotGrob(fig_a)
legend <- grob$grobs[[which(sapply(grob$grobs, function(x) x$name) == "guide-box")]]

# Combine figure A and B
combined_plot <- plot_grid(fig_a_nolegend, fig_b_nolegend, ncol = 2, labels = c("A", "B"))

# Add legend below the combined plot
final_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.05))
ggsave("results/figures/figs2.compare_yoy.tiff", final_plot, width = 13, height = 15)


# 5. Assess suitability for seasonal immunisation using 6-month moving interval approach----
# Threshold adjustment: 75% → 80%; proportional adjustment for MIA threshold: 65% → 69%
# Rationale: Adjusted thresholds are used to assess differences with observations from Guo et al.

# 5.1 Compare observed and predicted MIA clusters across 21 provinces----
# Predicted MIA clusters
MIA_predict_province <- predict_province %>%
  # Calculate 6-month rolling sum of AAP by province-month
  group_by(Province) %>%
  summarise(pred.AAP = list(c(pred.AAP, pred.AAP[1:5])), 
            Month = list(c(Month, 13:17))) %>%
  unnest(c(pred.AAP, Month)) %>%
  group_by(Province) %>%
  mutate(RollSum_6m = RcppRoll::roll_suml(pred.AAP, 6, fill = NA)) %>% 
  filter(Month <= 12) %>%
  # Compute AAP proportion captured by the 6-month window, and retain the maximum
  mutate(RollProp_6m = RollSum_6m / sum(pred.AAP) * 100) %>%
  filter(RollProp_6m == max(RollProp_6m)) %>%
  # Classify provinces as suitable for seasonal immunisation if ≥69% of AAP is concentrated in 6 months
  mutate(cluster_69 = as.factor(ifelse(RollProp_6m >= (65/75 * 80), 
                                       as.character(Month), "0"))) %>%   # 等比例调节cut-off，结果转换为factor
  dplyr::select(-pred.AAP)
# Add cluster labels and source
MIA_predict_province <- MIA_predict_province %>%
  mutate(cluster_label = case_when(
    cluster_69 == "9" ~ "September",
    cluster_69 == "10" ~ "October",
    cluster_69 == "11" ~ "November",
    cluster_69 == "12" ~ "December",
    cluster_69 == "0" ~ "Not suitable")
  ) %>%
  rename(cluster = cluster_69) %>%
  mutate(Source = "Predicted")

# Observed MIA clusters
MIA_observe_province <- observe_province %>%
  group_by(Province) %>%
  summarise(AAP = list(c(AAP, AAP[1:5])), 
            Month = list(c(Month, 13:17))) %>% 
  unnest(c(AAP, Month)) %>%
  group_by(Province) %>%
  mutate(RollSum_6m = RcppRoll::roll_suml(AAP, 6, fill = NA)) %>%
  filter(Month <= 12) %>%
  mutate(RollProp_6m = RollSum_6m / sum(AAP) * 100) %>%
  filter(RollProp_6m == max(RollProp_6m)) %>%
  mutate(cluster_80 = as.factor(ifelse(round(RollProp_6m, 0) >= 80,
                                       as.character(Month), "0"))) %>%
  dplyr::select(-AAP)
# Add cluster labels and source
MIA_observe_province <- MIA_observe_province %>%
  mutate(cluster_label = case_when(
    cluster_80 == "9" ~ "September",
    cluster_80 == "10" ~ "October",
    cluster_80 == "11" ~ "November",
    cluster_80 == "0" ~ "Not suitable")
  ) %>%
  rename(cluster = cluster_80) %>%
  mutate(Source = "Observed")

# Keep only 21 comparable provinces
MIA_compare <- bind_rows(
  MIA_observe_province, 
  MIA_predict_province %>% semi_join(MIA_observe_province, by = "Province")
  ) %>%
  arrange(Province)
write.csv(MIA_compare, "results/tables/tab1.MIA_cpmpare.csv")

# 5.2 Observed plus predicted MIA clusters across 31 provinces----
MIA_combined_province <- MIA_predict_province %>%
  anti_join(MIA_observe_province, by = c("Province")) %>%  # Only provinces without observed data
  rbind(MIA_observe_province)
MIA_combined_province_map <- left_join(MIA_combined_province, chinamap_ADM1, by = c("Province" = "name"))

# Filter provinces potentially suitable for seasonal immunisation: 23
MIA_combined_province_map %>%
  filter(cluster != "0") %>%
  pull(Province)

# Filter provinces with optimal beginning month of seasonal programme in October or November: 21
MIA_combined_province_map %>%
  filter(cluster %in% c("10","11")) %>%
  pull(Province)

# Figure 5
main_map2 <- ggplot() +

  # Base China provincial map
  geom_sf(data = chinamap_ADM1, color = "grey20", fill = "white", size = 5) +
  
  # First layer: fill colors for cluster
  geom_sf(data = st_as_sf(MIA_combined_province_map), 
          aes(fill = cluster_label), color = "grey20", alpha = 0.8) + 
  
  # Second layer: hatch patterns for source
  ggpattern::geom_sf_pattern(data = st_as_sf(MIA_combined_province_map), 
                  aes(pattern = Source),  # **只使用 pattern**
                  color = "grey20",  alpha = 0, 
                  pattern_density = 0.1, pattern_color = "black", 
                  pattern_fill = NA, pattern_spacing = 0.01, pattern_alpha = 0.5) +
  
  # Color legend
  scale_fill_manual(values = c("September" = "#508AB2", "October" = "#F6C14D", 
                               "November" = "#8768a6", "December" = "#9bbf8a", 
                               "Not suitable" = "#C52A20"),
                    breaks = c("September", "October", "November", "December", "Not suitable"),  
                    name = "Optimal beginning of \nseasonal immunisation") +   
  
  # Pattern legend
  scale_pattern_manual(values = c("Observed" = "none", 
                                  "Predicted" = "stripe"),
                       breaks = c("Observed", "Predicted"),
                       name = "Data source") +
  
  # Legend order
  guides(fill = guide_legend(order = 1),
         pattern = guide_legend(order = 2)) +
  
  # Reference latitude lines
  geom_hline(yintercept = 23.5, lty = 3, color = "darkgrey") +
  geom_hline(yintercept = 35, lty = 2, color = "darkgrey") +
  geom_text(aes(x = Inf, y = 23.5, label = "23.5°N"), hjust = 1, vjust = -1, size = 3, color = "grey40") +
  geom_text(aes(x = Inf, y = 35, label = "35°N"), hjust = 1, vjust = -1, size = 3, color = "grey40") +
  
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold")) +
  
  # Add red box around Hong Kong and Macao
  geom_rect(aes(xmin = hk_bbox["xmin"]-0.3, xmax = hk_bbox["xmax"]+0.1, 
                ymin = hk_bbox["ymin"]-0.3, ymax = hk_bbox["ymax"]+0.1), 
            color = adjustcolor("red2", alpha.f = 1), fill = NA, linewidth = 0.5)

# Macao inset map
macao_map2 <- ggplot() +
  geom_sf(data = st_as_sf(MIA_combined_province_map)%>% filter(Province == "Macao"), aes(fill = cluster_label), color = "grey20", alpha = 0.8) + 
  geom_sf_pattern(data = st_as_sf(MIA_combined_province_map) %>% filter(Province == "Macao"), aes(pattern = Source), 
                  color = "grey20", size = 5, alpha = 0,
                  pattern_density = 0.1, pattern_color = "black", 
                  pattern_fill = NA, pattern_spacing = 0.2, pattern_alpha = 0.5) +
  scale_fill_manual(values = c("September" = "#508AB2", "October" = "#F6C14D", 
                               "November" = "#8768a6", "December" = "#9bbf8a", 
                               "Not suitable" = "#C52A20")) +
  scale_pattern_manual(values = c("Observed" = "none", 
                                  "Predicted" = "stripe")) +
  theme_void() +
  theme(legend.position = "none")

# Hong Kong inset map
hk_map2 <- ggplot() +
  geom_sf(data = st_as_sf(MIA_combined_province_map)%>% filter(Province == "Hong Kong"), aes(fill = cluster_label), color = "grey20", alpha = 0.8) + 
  geom_sf_pattern(data = st_as_sf(MIA_combined_province_map) %>% filter(Province == "Hong Kong"), aes(pattern = Source), 
                  color = "grey20", size = 5, alpha = 0,
                  pattern_density = 0.1, pattern_color = "black", 
                  pattern_fill = NA, pattern_spacing = 0.2, pattern_alpha = 0.5) +
  scale_fill_manual(values = c("September" = "#508AB2", "October" = "#F6C14D", 
                               "November" = "#8768a6", "December" = "#9bbf8a", 
                               "Not suitable" = "#C52A20")) +
  scale_pattern_manual(values = c("Observed" = "none", 
                                  "Predicted" = "stripe")) + 
  theme_void() +
  theme(legend.position = "none")

# Combine main map and inset maps
final_map2 <- cowplot::ggdraw() +
  # Main map
  draw_plot(main_map2) +
  
  # Add inset boxes
  geom_rect(aes(xmin = 0.10, xmax = 0.36, ymin = 0.08, ymax = 0.32), 
            color = "red2", fill = "transparent", linewidth = 0.5) +  # 先画框
  
  # Insert Macao and Hong Kong maps
  draw_plot(macao_map2, x = 0.10, y = 0.15, width = 0.1, height = 0.1) +  
  draw_plot(hk_map2, x = 0.20, y = 0.15, width = 0.15, height = 0.15) + 
  
  # Connection lines
  geom_segment(aes(x = 0.458, y = 0.394, xend = 0.10, yend = 0.32),
               linetype = "dashed", color = "red2", alpha = 1, linewidth = 0.5) +  
  geom_segment(aes(x = 0.467, y = 0.38, xend = 0.36, yend = 0.08),
               linetype = "dashed", color = "red2", alpha = 1, linewidth = 0.5) +  
  
  # Labels for insets
  draw_text("Macao", x = 0.15, y = 0.12, size = 10, color = "black") +  
  draw_text("Hong Kong", x = 0.28, y = 0.12, size = 10, color = "black")  

# Save as TIFF (PDF causes hatch pattern issues)
ggsave("results/figures/fig5.combine_cluster.tiff", plot = final_map2, width = 7, height = 4,  bg = "white")
