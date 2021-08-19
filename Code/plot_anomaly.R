library(dplyr)
library(R.matlab)
library(sf)
library(plot.matrix)
library(colorspace)
library(optimbase)

# --- Source data binding function ---
source("./calanus_data/Code/bind_years.R")

# --- Load critical habitats (CH) ---
ch <- readMat("./calanus-for-whales/Code/CalWhaleRegions.mat")

# Roseway Basin
RB <- st_polygon(list(ch$Regions[,,2]$lonlat))
RB <- st_sf(st_sfc(RB))
# Cape Cod Bay
CCB <- st_polygon(list(ch$Regions[,,3]$lonlat))
CCB <- st_sf(st_sfc(CCB))
# Bay of Fundy
BOF <- st_polygon(list(ch$Regions[,,4]$lonlat))
BOF <- st_sf(st_sfc(BOF))
# Wilkinson Basin
WB <- st_polygon(list(ch$Regions[,,5]$lonlat))
WB <- st_sf(st_sfc(WB))
# Jordan Basin
JB <- st_polygon(list(ch$Regions[,,6]$lonlat))
JB <- st_sf(st_sfc(JB))
# Great South Channel
GSC <- st_polygon(list(ch$Regions[,,7]$lonlat))
GSC <- st_sf(st_sfc(GSC))
# Massachusetts Bay
MB <- st_polygon(list(ch$Regions[,,8]$lonlat))
MB <- st_sf(st_sfc(MB))

# ---- Load model data (md) ----
if (length(years) == 1) {
  md <- readr::read_csv(file.path(fp_md, paste0(years[1], ".csv"))) %>% 
    dplyr::filter(dataset == "ECOMON")
} else {
  md <- bind_years(fp = file.path(fp_md), years = years) %>%
    dplyr::filter(dataset == "ECOMON")
}

# -------- Take the log of count data --------
if (species == "cfin") {
  if (md$dataset == "ECOMON_STAGED") {
    md$abund <- as.data.frame(md[paste0(species, "_CV_VI")])$cfin_CV_VI + as.data.frame(md[paste0(species, "_CIV")])$cfin_CIV
  } else if (md$dataset == "ECOMON") {
    md$abund <- as.data.frame(md[paste0(species, "_CV_VI")])$cfin_CV_VI
  }
} else if (species == "ctyp") {
  md$abund <- as.data.frame(md[paste0(species, "_total")])$ctyp_total
} else if (species == "pcal") {
  md$abund <- as.data.frame(md[paste0(species, "_total")])$pcal_total
}


# --- Convert to sf to assign regions ---
md_sf <- st_as_sf(md, coords = c("lon", "lat"))

# --- Intersect with critical habitat polygons ---
iRB <- st_intersects(RB, md_sf, sparse = FALSE)
iCCB <- st_intersects(CCB, md_sf, sparse = FALSE)
iBOF <- st_intersects(BOF, md_sf, sparse = FALSE)
iWB <- st_intersects(WB, md_sf, sparse = FALSE)
iJB <- st_intersects(JB, md_sf, sparse = FALSE)
iGSC <- st_intersects(GSC, md_sf, sparse = FALSE)
iMB <- st_intersects(MB, md_sf, sparse = FALSE)

# --- Add names ---
poly_names <- rep("", nrow(md))
poly_names[iMB] <- "MB"
poly_names[iCCB] <- "CCB"
poly_names[iWB] <- "WB"
poly_names[iJB] <- "JB"
poly_names[iGSC] <- "GSC"

# --- Assign region names ---
md <- md %>%
  dplyr::mutate(region = poly_names)

# --- Assign NA to locations with no region ---
md$region[md$region == ""] <- NA

# --- Plot model data for Great South Channel ---

md_gsc <- md %>% dplyr::filter(region == "GSC") %>%
  dplyr::mutate(date = as.Date(paste(year, month, "01", sep = "-")))

md_gsc_clim <- md_gsc %>%
  dplyr::group_by(month) %>%
  summarize(climatological_mean = mean(abund, na.rm = TRUE))

ggplot(md_gsc, aes(x = date, y = abund)) +
  geom_point(alpha = 0.6) +
  ylim(c(0,6))

ggplot(md_gsc_clim, aes(x = as.factor(month), y = climatological_mean)) +
  geom_point()

thresholds <- c(10000)

md$thresh_10 <- if_else(md$abund > thresholds[1], 1, 0)

md_clim <- md %>%
  dplyr::group_by(month, region) %>%
  dplyr::summarize(climatological_mean = mean(abund, na.rm = TRUE),
                   threshold_mean = mean(thresh_10, na.rm = TRUE))

for (i in 1:12) {
  for (j in regions) {
      md_ij <- md %>% dplyr::filter(month == i,
                                    region == j,
                                    !is.na(region))
      md_clim_ij <- md_clim %>% dplyr::filter(month == i, 
                                              region == j,
                                              !is.na(region))
      if (!(NA %in% md_ij$abund)) {
        md$anomaly[md$month == i & md$region == j & !is.na(md$region)] <- md_ij$abund - md_clim_ij$climatological_mean
        md$thresh_anomaly[md$month == i & md$region == j & !is.na(md$region)] <- md_ij$thresh_10 - md_clim_ij$threshold_mean
      }
   }
}

# --- Initialize region names ---
regions <- c("CCB", "GSC", "MB", "JB", "WB")

# --- Compute anomaly and mean ---
md <- md %>% 
  group_by(year, region) %>%
  dplyr::summarize(mean_anomaly = mean(anomaly, na.rm = TRUE),
                   mean_thresh_anomaly = mean(thresh_anomaly, na.rm = TRUE))

# ---- 1,000 ----
# --- Initialize matrix ----
mat <- matrix(nrow = length(regions), ncol = length(years))
mat_thresh <- matrix(nrow = length(regions), ncol = length(years))

# --- Assign matrix values ----
for (k in 1:length(regions)) {
  for (l in 1:length(years)) {
    year_ <- years[l]
    region_ <- regions[k]
    year_md <- md %>% filter(year == year_, region == region_)
    mat[k,l] <- mean(year_md$mean_anomaly, na.rm = TRUE)
    mat_thresh[k,l] <- mean(year_md$mean_thresh_anomaly, na.rm = TRUE)
  }
}

# --- Remove NAs ---
mat[is.nan(mat)] <- NA
mat_thresh[is.nan(mat_thresh)] <- NA

temp_mat <- sign(mat)

mat <- log10(abs(mat)) * temp_mat

# --- Add row and column names ---
colnames(mat) <- years
rownames(mat) <- regions

colnames(mat_thresh) <- years
rownames(mat_thresh) <- regions

# --- Plot ---
par(mar=c(5, 5, 5, 10), cex.axis = 1.5, las = 1) 
plot(mat, xlab = "", ylab = "",
     fmt.key="%.1f",
     fmt.cell='%.1f',
     col = colorspace::diverge_hsv(11),
     main = "EcoMon Abundance Anomaly",
     digits = TRUE,
     text.cell = list(cex = 0.8))

par(mar=c(5, 5, 5, 10), cex.axis = 1.5, las = 1) 
plot(mat_thresh, xlab = "", ylab = "",
     fmt.key="%.1f",
     fmt.cell='%.1f',
     col = colorspace::diverge_hsv(11),
     main = "EcoMon Threshold Anomaly",
     digits = TRUE,
     text.cell = list(cex = 0.8))


corr_dat <- as.data.frame(pivot_longer(as.data.frame(mat), cols = 1:18))
corr_dat[is.na(corr_dat)] <- 0
corr_dat$name <- as.numeric(corr_dat$name)
corr_dat$region <- rep(c("CCB", "GSC", "MB", "JB", "WB"), 18)

corr_dat <- corr_dat %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(cor = stats::cor.test(name, value)$estimate,
                   sig = stats::cor.test(name, value)$p.value)

corr_dat <- as.data.frame(pivot_longer(as.data.frame(mat_thresh), cols = 1:18))
corr_dat[is.na(corr_dat)] <- 0
corr_dat$name <- as.numeric(corr_dat$name)
corr_dat$region <- rep(c("CCB", "GSC", "MB", "JB", "WB"), 18)

corr_dat <- corr_dat %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(cor = stats::cor.test(name, value)$estimate,
                   sig = stats::cor.test(name, value)$p.value)

# ---- 4,000 ----
# --- Initialize matrix ---
mat <- matrix(nrow = length(years), ncol = length(regions))

# --- Assign matrix values ----
for (k in 1:length(years)) {
  for (l in 1:length(regions)) {
    year_ <- years[k]
    region_ <- regions[l]
    year_md <- md %>% filter(year == year_, region == region_)
    mat[k,l] <- mean(year_md$mean_anomaly_thresh_4, na.rm = TRUE)
  }
}

# --- Remove NAs ---
mat[is.nan(mat)] <- NA
# --- Add row and column names ---
rownames(mat) <- years
colnames(mat) <- regions

# --- Plot ---
par(mar=c(5, 5, 5, 10), cex.axis = 1.5, las = 1) 
plot(mat, xlab = "", ylab = "",
     fmt.key="%.2e",
     col=colorspace::diverge_hsv(41),
     breaks = seq(-2/1e17, 2/1e17, by = 1e-18),
     main = expression(paste("4,000 m"^"-2")))

# ---- 10,000 ----
# --- Initialize matrix ---
mat <- matrix(nrow = length(years), ncol = length(regions))

# --- Assign matrix values ----
for (k in 1:length(years)) {
  for (l in 1:length(regions)) {
    year_ <- years[k]
    region_ <- regions[l]
    year_md <- md %>% filter(year == year_, region == region_)
    mat[k,l] <- mean(year_md$mean_anomaly_thresh_10, na.rm = TRUE)
  }
}

# --- Remove NAs ---
mat[is.nan(mat)] <- NA
# --- Add row and column names ---
rownames(mat) <- years
colnames(mat) <- regions

# --- Plot ---
par(mar=c(5, 5, 5, 10), cex.axis = 1.5, las = 1) 
plot(mat, xlab = "", ylab = "",
     fmt.key="%.2e",
     col=colorspace::diverge_hsv(41),
     breaks = seq(-2/1e17, 2/1e17, by = 1e-18),
     main = expression(paste("10,000 m"^"-2")))





