# Camille Ross
# 10 November 2020
# Purpose: Create table of zooplankton abundance per month per year

# Load libraries
library(readr)
library(dplyr)
library(ggplot2)

# Set working directory
DIR <- "~/Desktop/Calanus_Project/calanus_data/Data/Databases"
setwd(DIR)

# -------- EcoMon --------
# Load EcoMon zooplankton data
zoo_dat <- read_csv("zooplankton_database.csv") %>% 
  dplyr::filter(dataset == "ECOMON" & year >= 2000) %>%
  dplyr::select(month, year, cfin_CV_VI, ctyp_total, pcal_total)

# ---- Cfin ----
# Plot calanus finmarchicus late stage totals per month
ggplot(data = zoo_dat, mapping = aes(as.factor(month), cfin_CV_VI)) +
  stat_summary(fun = sum,
               geom = "bar")

# Plot calanus finmarchicus late stage totals per year
ggplot(data = zoo_dat, mapping = aes(as.factor(year), cfin_CV_VI)) +
  stat_summary(fun = sum,
               geom = "bar")

# ---- Ctyp ----
# Plot calanus finmarchicus late stage totals per month
ggplot(data = zoo_dat, mapping = aes(as.factor(month), ctyp_total)) +
  stat_summary(fun = sum,
               geom = "bar")

# Plot calanus finmarchicus late stage totals per year
ggplot(data = zoo_dat, mapping = aes(as.factor(year), ctyp_total)) +
  stat_summary(fun = sum,
               geom = "bar")

# ---- Pcal ----
# Plot calanus finmarchicus late stage totals per month
ggplot(data = zoo_dat, mapping = aes(as.factor(month), pcal_total)) +
  stat_summary(fun = sum,
               geom = "bar")

# Plot calanus finmarchicus late stage totals per year
ggplot(data = zoo_dat, mapping = aes(as.factor(year), pcal_total)) +
  stat_summary(fun = sum,
               geom = "bar")

# -------- CPR --------
# Load CPR zooplankton data
zoo_dat <- read_csv("zooplankton_database.csv") %>% 
  dplyr::filter(dataset == "CPR" & year >= 2000) %>%
  dplyr::select(month, year, cfin_CV_VI, ctyp_total, pcal_total)

# ---- Cfin ----
# Plot calanus finmarchicus late stage totals per month
ggplot(data = zoo_dat, mapping = aes(as.factor(month), cfin_CV_VI)) +
  stat_summary(fun = sum,
               geom = "bar")

# Plot calanus finmarchicus late stage totals per year
ggplot(data = zoo_dat, mapping = aes(as.factor(year), cfin_CV_VI)) +
  stat_summary(fun = sum,
               geom = "bar")

# ---- Ctyp ----
# Plot calanus finmarchicus late stage totals per month
ggplot(data = zoo_dat, mapping = aes(as.factor(month), ctyp_total)) +
  stat_summary(fun = sum,
               geom = "bar")

# Plot calanus finmarchicus late stage totals per year
ggplot(data = zoo_dat, mapping = aes(as.factor(year), ctyp_total)) +
  stat_summary(fun = sum,
               geom = "bar")

# ---- Pcal ----
# Plot calanus finmarchicus late stage totals per month
ggplot(data = zoo_dat, mapping = aes(as.factor(month), pcal_total)) +
  stat_summary(fun = sum,
               geom = "bar")

# Plot calanus finmarchicus late stage totals per year
ggplot(data = zoo_dat, mapping = aes(as.factor(year), pcal_total)) +
  stat_summary(fun = sum,
               geom = "bar")
  



