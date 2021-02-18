if (length(years) == 1) {
  md <- readr::read_csv(file.path(fp_md, paste0(years[1], ".csv"))) %>% 
    dplyr::filter(dataset %in% biomod_dataset)
} else {
  md <- bind_years(fp = file.path(fp_md), years = years) %>%
    dplyr::filter(dataset %in% biomod_dataset)
}

# -------- Compute anomaly --------
if (species == "cfin") {
  md <- md %>% dplyr::group_by(dataset) %>%
    dplyr::mutate(mean = mean(log10(`cfin_CV_VI` + 1), na.rm = TRUE),
                  sd = sd(log10(`cfin_CV_VI` + 1), na.rm = TRUE),
                  anomaly = (log10(`cfin_CV_VI` + 1) - mean) / sd) %>%
    dplyr::ungroup()
} else if (species == "ctyp") {
  md <- md %>% dplyr::group_by(dataset) %>%
    dplyr::mutate(mean = mean(log10(`ctyp_total` + 1), na.rm = TRUE),
                  sd = sd(log10(`ctyp_total` + 1), na.rm = TRUE),
                  anomaly = (log10(`ctyp_total` + 1) - mean) / sd) %>%
    dplyr::ungroup()
} else if (species == "pcal") {
  md <- md %>% dplyr::group_by(dataset) %>%
    dplyr::mutate(mean = mean(log10(`pcal_total` + 1), na.rm = TRUE),
                  sd = sd(log10(`pcal_total` + 1), na.rm = TRUE),
                  anomaly = (log10(`pcal_total` + 1) - mean) / sd) %>%
    dplyr::ungroup()
}

# -------- Take the log of count data --------
if (species == "cfin") {
  md$abund <- as.data.frame(md[paste0(species, "_CV_VI")])$cfin_CV_VI
} else if (species == "ctyp") {
  md$abund <- as.data.frame(md[paste0(species, "_total")])$ctyp_total
} else if (species == "pcal") {
  md$abund <- as.data.frame(md[paste0(species, "_total")])$pcal_total
}

# -------- Exclude NAs and select columns --------
md <- md %>% dplyr::select(lat, lon, year, month, abund, wind, fetch, jday, chl, int_chl, bots, bott, sss, sst, sst_grad, lag_sst, uv, bat, dist, slope) %>%
  as.data.frame() %>%
  na.exclude() %>%
  dplyr::mutate(season = if_else(month %in% c(1:3), 1,
                                 if_else(month %in% c(4:6), 2,
                                         if_else(month %in% c(7:9), 3, 4))),
                region = if_else(lat <= 41.5 & lon < -70, "MAB", 
                                 if_else(lat >= 39 & lat <= 42 & lon >= -70 & lon <= -68, "GBK", "GOM")))


mab_mat <- matrix(nrow = length(years), ncol = 12)

mab <- md %>% dplyr::filter(region == "MAB")

for (i in 1:length(years)) {
  for (j in 1:12) {
    year_ <- years[i]
    month_md <- mab %>% filter(year == year_, month == j)
    month_md$abund <- if_else(month_md$abund > threshold, 1, 0)
    mab_mat[i,j] <- mean(month_md$abund)
  }
}

mab_mat[is.nan(mab_mat)] <- NA
rownames(mab_mat) <- years

par(mar=c(5.1, 4.1, 4.1, 4.1)) 
plot(mab_mat, xlab = "Month", ylab = "Year", main = "Mid-Atlantic Bight", key = NULL,
     col=inferno(10), breaks=c(0, 0.1, 0.2, 0.3, 0.4, 
                               0.5, 0.6, 0.7, 0.8, 0.9, 1))


gbk_mat <- matrix(nrow = length(years), ncol = 12)

gbk <- md %>% dplyr::filter(region == "GBK")

for (i in 1:length(years)) {
  for (j in 1:12) {
    year_ <- years[i]
    month_md <- gbk %>% filter(year == year_, month == j)
    month_md$abund <- if_else(month_md$abund > threshold, 1, 0)
    gbk_mat[i,j] <- mean(month_md$abund)
  }
}

gbk_mat[is.nan(gbk_mat)] <- NA
rownames(gbk_mat) <- years

par(mar=c(5.1, 4.1, 4.1, 4.1)) 
plot(gbk_mat, xlab = "Month", ylab = "Year", main = "George's Bank", key= NULL,
     col=inferno(10), breaks=c(0, 0.1, 0.2, 0.3, 0.4, 
                               0.5, 0.6, 0.7, 0.8, 0.9, 1))


gom_mat <- matrix(nrow = length(years), ncol = 12)

gom <- md %>% dplyr::filter(region == "GOM")

for (i in 1:length(years)) {
  for (j in 1:12) {
    year_ <- years[i]
    month_md <- gom %>% filter(year == year_, month == j)
    month_md$abund <- if_else(month_md$abund > threshold, 1, 0)
    gom_mat[i,j] <- mean(month_md$abund)
  }
}

gom_mat[is.nan(gom_mat)] <- NA
rownames(gom_mat) <- years

par(mar=c(5.1, 4.1, 4.1, 4.1)) 
plot(gom_mat, xlab = "Month", ylab = "Year", main = "Gulf of Maine", key = NULL,
     col=inferno(10), breaks=c(0, 0.1, 0.2, 0.3, 0.4, 
                               0.5, 0.6, 0.7, 0.8, 0.9, 1))





