#!/usr/bin/env Rscript

# R script
# Environmental conditions before, during and after the Orange County oil spill, 2021.

# This script contains R code for plotting environmental data, including mean wave direction,
# surface water temperature, salinity, chlorophyll a, CUTI, BEUTI, PIC and PON

# R version: 4.3.0

library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(tidyr)
library(ggpubr)




# ================================================================================================================ #
#### 1. NOAA 46253 ####

# Data for Mean Wave Direction

##### 1.1. Import data #####

# Define column names
cols <- c("YY","MM","DD","hh","mm","WDIR","WSPD","GST","WVHT","DPD","APD",
          "MWD","PRES","ATMP","WTMP","DEWP","VIS","TIDE")

# Read file
noaa_46253 <- read_table("noaa_46253_mwd.txt",
  skip = 2,
  col_names = cols,
  col_types = cols(
    YY = col_integer(),
    MM = col_integer(),
    DD = col_integer(),
    hh = col_integer(),
    mm = col_integer(),
    WDIR = col_double(),
    WSPD = col_double(),
    GST  = col_double(),
    WVHT = col_double(),
    DPD  = col_double(),
    APD  = col_double(),
    MWD  = col_double(),
    PRES = col_double(),
    ATMP = col_double(),
    WTMP = col_double(),
    DEWP = col_double(),
    VIS  = col_double(),
    TIDE = col_double()))


# Helpers for figures

# Format dates
noaa_46253 <- noaa_46253 %>%
  mutate(date = make_date(year = YY, month = MM, day = DD))

# Key sampling dates for figures (time)
time_days <- tibble(date = ymd(c("2021-10-19", "2021-11-04", "2021-11-19", "2021-12-06")),
                  label = c("18d", "34d", "49d", "66d"))
oil_spill <- ymd("2021-10-01")

# Vertical bars to identify each time
bars_time <- time_days %>%
  mutate(label = factor(label, levels = c("18d","34d","49d","66d"))) %>%
  transmute(date, xint = as.numeric(date), label)

# Vertical bars to identify date of oil spill
bar_oil  <- tibble(xint = as.numeric(oil_spill)) 
label_oil  <- tibble(date = oil_spill, xint = as.numeric(oil_spill), label = "Spill")

# Color palette
oil_palette <- colorRampPalette(c("gray40", rev(brewer.pal(4, "Blues")), "#1dd3b0"))(4)
oil_palette_bars <- setNames(oil_palette, levels(bars_time$label)) # color for vertical bars


##### 1.2. Mean Wave Direction ####

# Filter Oct–Dec and summarize
oct_dec_mwd <- noaa_46253 %>%
  filter(date >= ymd("2021-09-15"), date <= ymd("2021-12-15")) %>%
  group_by(date) %>%
  summarize(MWD_daily = mean(MWD, na.rm = TRUE), .groups = "drop")

# Plot
ymin_lim <- 150
ymax_lim <- 300

mwd_46253 <- ggplot(oct_dec_mwd, aes(x = date, y = MWD_daily)) +
  geom_vline(data = bars_time, aes(xintercept = xint, color = label),
             linetype = "dotted", linewidth = 0.3, inherit.aes = FALSE, show.legend = FALSE) +
  geom_vline(data = bar_oil, aes(xintercept = xint),
             linetype = "solid", color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  geom_line(color = "gray40", linetype = "dashed", linewidth = 0.3) +
  geom_point(color = "gray40", fill = "white", size = 0.7, shape = 21) +
  geom_text(data = bars_time, aes(x = date, label = as.character(label)),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_text(data = label_oil, aes(x = date, label = label),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_point(data = time_days %>% left_join(oct_dec_mwd, by = "date"), 
             aes(x = date, y = MWD_daily, fill = label),
             shape = 21, color = "black", size = 2.5, inherit.aes = FALSE) +
  scale_fill_manual(values = oil_palette, guide = "none") +
  scale_color_manual(values = oil_palette_bars, guide = "none") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(limits = c(ymin_lim, ymax_lim), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Date", y = "Mean Wave Direction (°)") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        legend.position = "none",
        plot.margin = margin(10, 10, 20, 10))

plot(mwd_46253)




# ================================================================================================================ #
#### 2. JWS-HB-Buoy ####

##### 2.1. Import data #####

# Read file
jws <- read_csv("jws_hb_buoy_water_temp.csv",
                col_types = cols(
                  time = col_datetime(format = ""),
                  sea_water_temperature = col_double(),
                  z = col_double()))

# Extract date
jws <- jws %>%
  mutate(date = as_date(time))


##### 3.2. Surface Water Temperature #####

# Filter Oct–Dec and summarize
jws_daily_temp <- jws %>%
  filter(date >= ymd("2021-09-15"), date <= ymd("2021-12-15")) %>%
  group_by(date) %>%
  summarize(WT_daily = mean(sea_water_temperature, na.rm = TRUE), .groups = "drop")

# Plot
ymin_lim <- min(jws_daily_temp$WT_daily, na.rm = TRUE)
ymax_lim <- max(jws_daily_temp$WT_daily, na.rm = TRUE)

water_temp_jws <- ggplot(jws_daily_temp, aes(x = date, y = WT_daily)) +
  geom_vline(data = bars_time, aes(xintercept = xint, color = label),
             linetype = "dotted", linewidth = 0.3, inherit.aes = FALSE, show.legend = FALSE) +
  geom_vline(data = bar_oil, aes(xintercept = xint),
             linetype = "solid", color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  geom_line(color = "gray40", linetype = "dashed", linewidth = 0.3) +
  geom_point(color = "gray40", fill = "white", size = 0.7, shape = 21) +
  geom_text(data = bars_time, aes(x = date, label = as.character(label)),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_text(data = label_oil, aes(x = date, label = label),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_point(
    data = time_days %>% dplyr::left_join(jws_daily_temp, by = "date"),
    aes(x = date, y = WT_daily, fill = label),
    shape = 21, color = "black", size = 2.5, inherit.aes = FALSE) +
  scale_fill_manual(values = oil_palette, guide = "none") +
  scale_color_manual(values = oil_palette_bars, guide = "none") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(limits = c(ymin_lim, ymax_lim), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Date", y = "Surface Water Temp. (°C)") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 7), 
    legend.position = "none",
    plot.margin = margin(10, 10, 20, 10))

plot(water_temp_jws)





# ================================================================================================================ #
#### 4. SCCOOS ####


##### 4.1. Salinity #####

# Read file
sccoos_sal <- read_csv("sccoos_newport_salinity.csv",
                       col_types = cols(
                         time = col_datetime(format = ""),
                         sea_water_practical_salinity_ctd = col_double(),
                         sea_water_practical_salinity_ctd_qc_agg = col_integer(),
                         z = col_double())) %>%
  filter(sea_water_practical_salinity_ctd_qc_agg == 1) %>% # Keep only good quality (flag == 1)
  mutate(date = as_date(time))

# Filter Oct–Dec and summarize
sccoos_sal_daily <- sccoos_sal %>%
  filter(date >= ymd("2021-09-15"), date <= ymd("2021-12-15")) %>%
  group_by(date) %>%
  summarize(SAL_daily = mean(sea_water_practical_salinity_ctd, na.rm = TRUE), .groups = "drop")

# Create a complete daily sequence and insert NAs where data are missing
sccoos_sal_daily_full <- sccoos_sal_daily %>%
  complete(date = seq(min(date), max(date), by = "1 day"))

# Plot
ymin_lim <- min(sccoos_sal_daily_full$SAL_daily, na.rm = TRUE)
ymax_lim <- max(sccoos_sal_daily_full$SAL_daily, na.rm = TRUE)

salinity_sccoos <- ggplot(sccoos_sal_daily_full, aes(x = date, y = SAL_daily)) +
  geom_vline(data = bars_time, aes(xintercept = xint, color = label),
             linetype = "dotted", linewidth = 0.3, inherit.aes = FALSE, show.legend = FALSE) +
  geom_vline(data = bar_oil, aes(xintercept = xint),
             linetype = "solid", color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  geom_line(color = "gray40", linetype = "dashed", linewidth = 0.3, na.rm = FALSE) +
  geom_point(color = "gray40", fill = "white", size = 0.7, shape = 21, na.rm = TRUE) +
  geom_text(data = bars_time, aes(x = date, label = as.character(label)),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_text(data = label_oil, aes(x = date, label = label),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_point(data = time_days %>% dplyr::left_join(sccoos_sal_daily_full, by = "date"),
             aes(x = date, y = SAL_daily, fill = label),
             shape = 21, color = "black", size = 2.5, inherit.aes = FALSE) +
  scale_fill_manual(values = oil_palette, guide = "none") +
  scale_color_manual(values = oil_palette_bars, guide = "none") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(limits = c(ymin_lim, ymax_lim), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Date", y = "Surface Salinity (PSS)") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 7),
    legend.position = "none",
    plot.margin = margin(10, 10, 20, 10))

plot(salinity_sccoos)




##### 4.2. Chlorophyll #####

# Read file
sccoos_chl <- read_csv("sccoos_newport_chlorophyll.csv",
                na = c("NaN", "NA", ""),
                col_types = cols(
                  time = col_datetime(format = ""),
                  mass_concentration_of_chlorophyll_in_sea_water_ctd = col_double(),
                  mass_concentration_of_chlorophyll_in_sea_water_ctd_qc_agg = col_double(),
                  z = col_double())) %>%
  rename(chl = mass_concentration_of_chlorophyll_in_sea_water_ctd,
         qc  = mass_concentration_of_chlorophyll_in_sea_water_ctd_qc_agg) %>%
  filter(qc == 1) %>%                 # keep only good-quality data
  mutate(date = as_date(time))


# Filter Oct–Dec and summarize
sccoos_chl_daily <- sccoos_chl %>%
  filter(date >= ymd("2021-09-15"), date <= ymd("2021-12-15")) %>%
  group_by(date) %>%
  summarize(CHL_daily = mean(chl, na.rm = TRUE), .groups = "drop")


# Plot
ymin_lim <- min(sccoos_chl_daily$CHL_daily, na.rm = TRUE)
ymax_lim <- max(sccoos_chl_daily$CHL_daily, na.rm = TRUE)

chlorophyll_sccoos <- ggplot(sccoos_chl_daily, aes(x = date, y = CHL_daily)) +
  # sampling bars (colored by label)
  geom_vline(data = bars_time, aes(xintercept = xint, color = label),
             linetype = "dotted", linewidth = 0.3, inherit.aes = FALSE, show.legend = FALSE) +
  # oil spill bar (black)
  geom_vline(data = bar_oil, aes(xintercept = xint),
             linetype = "solid", color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  # time series (solid gray line + hollow points)
  geom_line(color = "gray40", linetype = "dashed", linewidth = 0.3) +
  geom_point(color = "gray40", fill = "white", size = 0.7, shape = 21) +
  # labels pinned to the top (won’t shift with y-limits)
  geom_text(data = bars_time, aes(x = date, label = as.character(label)),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_text(data = label_oil, aes(x = date, label = label),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  # sample dots (filled with your oil_palette, outlined black)
  geom_point(
    data = time_days %>% dplyr::left_join(sccoos_chl_daily, by = "date"),
    aes(x = date, y = CHL_daily, fill = label),
    shape = 21, color = "black", size = 2.5, inherit.aes = FALSE) +
  scale_fill_manual(values = oil_palette, guide = "none") +
  scale_color_manual(values = oil_palette_bars, guide = "none") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(limits = c(ymin_lim, ymax_lim), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Date", y = "Chlorophyll (µg/L)") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        legend.position = "none",
        plot.margin = margin(10, 10, 20, 10))

plot(chlorophyll_sccoos)





# ================================================================================================================ #
#### 5. UPWELLING ####

##### 5.1. CUTI ####

# Read file
cuti <- read_csv("CUTI_daily.csv",
                 col_types = cols(
                   year  = col_integer(),
                   month = col_integer(),
                   day   = col_integer(),
                   .default = col_double())) %>%
  mutate(date = make_date(year, month, day))

# Filter Oct–Dec
cuti_daily <- cuti %>%
  transmute(date, CUTI = `33N`) %>%                     
  filter(date >= ymd("2021-09-15"), date <= ymd("2021-12-15"))

# Plot
ymin_lim <- min(cuti_daily$CUTI, na.rm = TRUE)
ymax_lim <- max(cuti_daily$CUTI, na.rm = TRUE)

cuti_plot <- ggplot(cuti_daily, aes(x = date, y = CUTI)) +
  geom_vline(data = bars_time, aes(xintercept = xint, color = label),
             linetype = "dotted", linewidth = 0.3, inherit.aes = FALSE, show.legend = FALSE) +
  geom_vline(data = bar_oil, aes(xintercept = xint),
             linetype = "solid", color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  geom_line(color = "gray40", linetype = "dashed", linewidth = 0.3) +
  geom_point(color = "gray40", fill = "white", size = 0.7, shape = 21) +
  geom_text(data = bars_time, aes(x = date, label = as.character(label)),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_text(data = label_oil, aes(x = date, label = label),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  # sample dots (filled by palette, black outline)
  geom_point(data = time_days %>% dplyr::left_join(cuti_daily, by = "date"),
             aes(x = date, y = CUTI, fill = label),
             shape = 21, color = "black", size = 2.5, inherit.aes = FALSE) +
  scale_fill_manual(values = oil_palette, guide = "none") +
  scale_color_manual(values = oil_palette_bars, guide = "none") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(limits = c(ymin_lim, ymax_lim), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Date", y = "CUTI") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        legend.position = "none",
        plot.margin = margin(10, 10, 20, 10))

plot(cuti_plot)



##### 5.2. BEUTI ####

# Read file
beuti <- read_csv("BEUTI_daily.csv",
                  col_types = cols(
                    year  = col_integer(),
                    month = col_integer(),
                    day   = col_integer(),
                    .default = col_double())) %>%
  mutate(date = make_date(year, month, day))

# Filter Oct–Dec
beuti_sel <- beuti %>%
  transmute(date, BEUTI = `33N`) %>%      
  filter(date >= ymd("2021-09-15"), date <= ymd("2021-12-15"))

# Plot
ymin_lim <- min(beuti_sel$BEUTI, na.rm = TRUE)
ymax_lim <- max(beuti_sel$BEUTI, na.rm = TRUE)

beuti_plot <- ggplot(beuti_sel, aes(x = date, y = BEUTI)) +
   geom_vline(data = bars_time, aes(xintercept = xint, color = label),
             linetype = "dotted", linewidth = 0.3, inherit.aes = FALSE, show.legend = FALSE) +
   geom_vline(data = bar_oil, aes(xintercept = xint),
             linetype = "solid", color = "black", linewidth = 0.3, inherit.aes = FALSE) +
   geom_line(color = "gray40", linetype = "dashed", linewidth = 0.3) +
  geom_point(color = "gray40", fill = "white", size = 0.7, shape = 21) +
   geom_text(data = bars_time, aes(x = date, label = as.character(label)),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_text(data = label_oil, aes(x = date, label = label),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
   geom_point(data = time_days %>% dplyr::left_join(beuti_sel, by = "date"),
             aes(x = date, y = BEUTI, fill = label),
             shape = 21, color = "black", size = 2.5, inherit.aes = FALSE) +
  scale_fill_manual(values = oil_palette, guide = "none") +
  scale_color_manual(values = oil_palette_bars, guide = "none") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(limits = c(ymin_lim, ymax_lim), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Date", y = "BEUTI") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        legend.position = "none",
        plot.margin = margin(10, 10, 20, 10))

plot(beuti_plot)




# ================================================================================================================ #
#### 6. BCO-DMO ####

# Read file
bco <- read_csv(
  "bco_dmo_pom.csv",
  col_select = c(Date, POC_Avg, PON_Avg),
  guess_max = 1e6, show_col_types = FALSE) %>%
  mutate(Date_chr = as.character(Date),
         date = ymd(Date_chr, quiet = TRUE),
         date = if_else(is.na(date), mdy(Date_chr, quiet = TRUE), date),
         date = if_else(is.na(date) & grepl("^[0-9]+$", Date_chr),
                   as.Date(as.numeric(Date_chr), origin = "1899-12-30"),
                   date)) %>% arrange(date)

# Filter Oct–Dec
bco_pom <- bco %>%
  filter(date >= ymd("2021-09-15"), date <= ymd("2021-12-15"))


##### 6.4. POC_Avg #####

# Plot 
ymin_lim <- min(bco_pom$POC_Avg, na.rm = TRUE)
ymax_lim <- max(bco_pom$POC_Avg, na.rm = TRUE)

bco_poc <- ggplot(bco_pom, aes(x = date, y = POC_Avg)) +
  geom_vline(data = bars_time, aes(xintercept = xint, color = label),
             linetype = "dotted", linewidth = 0.3, inherit.aes = FALSE, show.legend = FALSE) +
  geom_vline(data = bar_oil, aes(xintercept = xint),
             linetype = "solid", color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  geom_line(data = dplyr::filter(bco_pom, !is.na(POC_Avg)),
    aes(group = 1), color = "gray40", linetype = "dashed", linewidth = 0.3) +
  geom_point(color = "gray40", fill = "white", size = 0.7, shape = 21, na.rm = TRUE) +
  geom_text(data = bars_time, aes(x = date, label = as.character(label)),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_text(data = label_oil, aes(x = date, label = label),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = oil_palette, guide = "none") +
  scale_color_manual(values = oil_palette_bars, guide = "none") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(limits = c(ymin_lim, ymax_lim), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Date", y = "POC (\u03BCM)") +
  coord_cartesian(clip = "off") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        legend.position = "none",
        plot.margin = margin(10, 10, 20, 10))

plot(bco_poc)




##### 6.5. PON_Avg #####

# Plot 
ymin_lim <- min(bco_pom$PON_Avg, na.rm = TRUE)
ymax_lim <- max(bco_pom$PON_Avg, na.rm = TRUE)

bco_pon <- ggplot(bco_pom, aes(x = date, y = PON_Avg)) +
  geom_vline(data = bars_time, aes(xintercept = xint, color = label),
             linetype = "dotted", linewidth = 0.3, inherit.aes = FALSE, show.legend = FALSE) +
  geom_vline(data = bar_oil, aes(xintercept = xint),
             linetype = "solid", color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  geom_line(data = dplyr::filter(bco_pom, !is.na(PON_Avg)),
            aes(group = 1), color = "gray40", linetype = "dashed", linewidth = 0.3) +
  geom_point(color = "gray40", fill = "white", size = 0.7, shape = 21, na.rm = TRUE) +
  geom_text(data = bars_time, aes(x = date, label = as.character(label)),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_text(data = label_oil, aes(x = date, label = label),
            y = Inf, vjust = 1.7, hjust = 1.2,
            color = "black", fontface = "bold", size = 2.5,
            inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = oil_palette, guide = "none") +
  scale_color_manual(values = oil_palette_bars, guide = "none") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(limits = c(ymin_lim, ymax_lim), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Date", y = "PON (\u03BCM)") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        legend.position = "none",
        plot.margin = margin(10, 10, 20, 10))

plot(bco_pon)




#### 7. FIGURE ARTICLE ####

##### Figure 4 #####

# Modify border
lock_border <- theme(
  axis.title = element_text(size = 9, face = "bold"),
  panel.border = element_rect(color = "grey55", fill = NA, linewidth = 0.5),
  panel.background = element_rect(fill = "white", color = NA))

cuti_plotB          <- cuti_plot + lock_border
beuti_plotB         <- beuti_plot + lock_border
water_temp_jwsB     <- water_temp_jws     + scale_y_continuous(limits = c(0, NA)) + lock_border
bco_pocB            <- bco_poc            + scale_y_continuous(limits = c(0, NA)) + lock_border
salinity_sccoosB    <- salinity_sccoos    + scale_y_continuous(limits = c(33, NA)) + lock_border
bco_ponB            <- bco_pon            + scale_y_continuous(limits = c(0, NA)) + lock_border
mwd_46253B          <- mwd_46253          + scale_y_continuous(limits = c(100, NA)) + lock_border
chlorophyll_sccoosB <- chlorophyll_sccoos + scale_y_continuous(limits = c(0, NA)) + lock_border


# Create plot
ggarrange(cuti_plotB,beuti_plotB,
          water_temp_jwsB, bco_pocB,
          salinity_sccoosB, bco_ponB,
          mwd_46253B, chlorophyll_sccoosB, 
          ncol = 2, nrow = 4,                  
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"), 
          font.label = list(size = 10),
          common.legend = T,               
          legend = "top")

ggsave("Figure5.tiff", last_plot(), width = 7, height = 9, dpi = 600)

