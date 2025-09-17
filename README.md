# RSV Seasonality Modelling in China
 
This repository contains the data and core codes for the paper: **Leveraging meteorological data for understanding regional distribution of RSV seasonality and its implication for infant immunisation strategy in China: a modelling study**

**Authors:** Shuyu Deng, Ling Guo, Shiqi Sun, Xin Wang, You Li.

**Corresponding author:**
Prof. You Li (You.Li@njmu.edu.cn), Prof. Xin Wang (Xin.Wang@njmu.edu.cn).

---

## Overview

All data processing, analysis, and output generation can be run from 'core.R'. This script loads all dependencies and custom functions from other scripts and executes the full workflow. The repository is structured to allow reproducible analyses from raw RSV and meteorological data to final manuscript results.

---

## Repository structure

```
RSV_seasonality_China_R/
├── RSV_seasonality_China_R.Rproj  # RStudio project file
├── core.R                         # Main script to run full analysis and generate results
├── dataProcessing.R               # RSV and meteorological data cleaning and processing
├── descriptiveFun.R               # Custom functions for descriptive statistics
├── weatherFun.R                   # Custom functions for meteorological data processing
├── modelFun.R                     # Custom functions for model fitting, prediction, and plotting
├── data/                          # Input datasets
│   ├── chinamap/                  # Shapefiles or JSON for basic China maps
│   │   ├── ADM1/                  # National-level shapefiles
│   │   └── ADM2/                  # Province-level JSON
│   ├── rsvdata/                   # RSV activity datasets
│   └── weatherdata/               # Meteorological datasets
├── results/                       # Output results
│   ├── figures/                   # Figures for manuscript
│   └── tables/                    # Tables and other results for manuscript
```

---

## System requirements

* **Operating system:** macOS Sonoma 14.1.1 / Windows 10 or later (tested on MacBook Pro 14", 2023, Apple M2 Pro, 16GB RAM)
* **Software:**
  * R version 4.3.2
  * RStudio (any recent version)
* **R packages used with versions:**

```
imputeTS_3.3   geosphere_1.5-20   GSODR_4.1.3   wktmo_1.0.5   lubridate_1.9.4   aweek_1.0.3
purrr_1.0.4    tidyr_1.3.1       dplyr_1.1.4   readxl_1.4.5  ggpattern_1.1.4   cowplot_1.1.3
RcppRoll_0.3.1 MASS_7.3-65       jsonlite_2.0.0 sf_1.0-21     ggplot2_3.5.2     stringr_1.5.1
```

* No non-standard hardware is required.

---

## Installation guide

1. Install R and RStudio.

2. Clone the repository:

```bash
git clone https://github.com/ShuyuDeng/RSV_seasonality_China_R.git
```

3. Open `RSV_seasonality_China_R.Rproj` in RStudio.

4. Install required R packages (if not already installed):

```r
install.packages(c("readxl","dplyr","tidyr","purrr","aweek","lubridate","wktmo","GSODR",
                   "geosphere","imputeTS","stringr","ggplot2","sf","jsonlite","MASS",
                   "RcppRoll","cowplot","ggpattern"))
```

5. Typical install time: ~5–10 minutes on a standard desktop/laptop with internet connection.

---

## Demo and instructions for use

* **Input data:** All raw RSV and meteorological data are provided in the `data/` folder.

* **Running demo:** Open `core.R` in RStudio and run it **line by line** to execute the full workflow. This includes RSV data cleaning, meteorological data processing, model fitting, prediction, and figure/table generation. **Do not run the script with `source()`**, as execution may pause partway through.

* **Expected output:** Results saved in the `results/` folder.

* **Expected run time:** Less than 5 minutes on a normal desktop computer. Since processing of meteorological data and some analyses can be time-consuming, processed files have already been stored in `data/weatherdata/` and `results/`. This substantially reduces the overall run time. To regenerate all intermediate data or results, set the custom function parameter `redo = TRUE` (default is FALSE).

* **Reproducing manuscript results:** The repository contains all raw data and core codes used in the study. Running `core.R` with the provided data will reproduce all results reported in the manuscript.

---

## License

This repository is licensed under the **GNU General Public License v3.0**.
