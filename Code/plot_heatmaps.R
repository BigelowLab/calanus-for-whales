library(plot.matrix)
library(R.matlab)
library(sf)
library(plyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(marmap)
library(ggnewscale)
library(tidyr)
library(raster)

# ---- Load critical habitats (CH) ----
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

if (length(years) == 1) {
  md <- readr::read_csv(file.path(fp_md, paste0(years[1], ".csv"))) %>% 
    dplyr::filter(dataset %in% "ECOMON")
} else {
  md <- bind_years(fp = file.path(fp_md), years = years) %>%
    dplyr::filter(dataset %in% "ECOMON")
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
} else if (species == "pseudo") {
  md <- md %>% dplyr::group_by(dataset) %>%
    dplyr::mutate(mean = mean(log10(`pseudo_total` + 1), na.rm = TRUE),
                  sd = sd(log10(`pseudo_total` + 1), na.rm = TRUE),
                  anomaly = (log10(`pseudo_total` + 1) - mean) / sd) %>%
    dplyr::ungroup()
}

# -------- Take the log of count data --------
if (species == "cfin") {
  md$abund <- as.data.frame(md[paste0(species, "_CV_VI")])$cfin_CV_VI
} else if (species == "ctyp") {
  md$abund <- as.data.frame(md[paste0(species, "_total")])$ctyp_total
} else if (species == "pseudo") {
  md$abund <- as.data.frame(md[paste0(species, "_total")])$pseudo_total
}

# -------- Exclude NAs and select columns --------
md <- md %>% 
  as.data.frame() %>%
  na.exclude() %>%
  dplyr::mutate(season = if_else(month %in% c(1:3), 1,
                                 if_else(month %in% c(4:6), 2,
                                         if_else(month %in% c(7:9), 3, 4))))


md_sf <- st_as_sf(md, coords = c("lon", "lat"))

iRB <- st_intersects(RB, md_sf, sparse = FALSE)
iCCB <- st_intersects(CCB, md_sf, sparse = FALSE)
iBOF <- st_intersects(BOF, md_sf, sparse = FALSE)
iWB <- st_intersects(WB, md_sf, sparse = FALSE)
iJB <- st_intersects(JB, md_sf, sparse = FALSE)
iGSC <- st_intersects(GSC, md_sf, sparse = FALSE)
iMB <- st_intersects(MB, md_sf, sparse = FALSE)

poly_names <- rep("", nrow(md))
poly_names[iMB] <- "MB"
poly_names[iCCB] <- "CCB"
poly_names[iWB] <- "WB"
poly_names[iJB] <- "JB"
poly_names[iGSC] <- "GSC"

md <- md %>%
  dplyr::mutate(region = poly_names)

md$region[md$region == ""] <- NA

blues <- c("lightsteelblue4", "lightsteelblue3",
           "lightsteelblue2", "lightsteelblue1")
bat <- getNOAA.bathy(lon1 = -76, lon2 = -64,
              lat1 = 35, lat2 = 47, resolution = 10)

plot(bat, image = TRUE, land = TRUE, lwd = 0.1,
     bpal = list(c(0, max(bat), "grey"),
                 c(min(bat),0,blues)))

ggplot(bat, aes(x=x, y=y)) + 
  # Set lat/lon bounds
  coord_quickmap() +
  # Background
  geom_raster(aes(fill=z)) +
  # Add fill colors
  scale_fill_etopo() +
  # Contours
  geom_contour(aes(z=z),
               breaks=c(0, -100, -200, -500, -1000, -2000, -4000),
               colour="black", size=0.2) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("") + ylab("") + 
  geom_point(data = md, aes(lon, lat), pch = 16, cex = 1) +
  new_scale_fill() +
  # Add projection data
  geom_tile(data = md %>% filter(!is.na(region)), aes(lon, lat, fill = region)) +
  # Add projection color gradient and label
  scale_fill_viridis_d(na.value = "white") +
  labs(x = "",
       y = "",
       fill = "Region") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank())

ggplot() +
  geom_point(data = md, aes(lon, lat), pch = 1, cex = 2) +
  # Add projection data
  geom_tile(data = md %>% filter(!is.na(region)), aes(lon, lat, fill = region)) +
  # Add projection color gradient and label
  scale_fill_viridis_d(na.value = "white") +
  labs(x = "",
       y = "") +
  # Add world map data
  geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
  coord_quickmap(xlim = c(-76, -64),
                 ylim = c(35, 47),
                 expand = TRUE) +
  # Remove grid lines
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_heatmap <- function(region_, region_name_, years, md, td) {
  
  mat <- matrix(nrow = length(years), ncol = 12)
  region_md <- md %>% dplyr::filter(region == region_)
  
  for (i in 1:length(years)) {
    for (j in 1:12) {
      year_ <- years[i]
      month_md <- region_md %>% filter(year == year_, month == j)
      month_md$abund <- if_else(month_md$abund > td, 1, 0)
      mat[i,j] <- mean(month_md$abund)
    }
  }
  
  mat[is.nan(mat)] <- NA
  rownames(mat) <- years
  
  par(mfrow = c(5, 3), mar=c(5.1, 4.1, 4.1, 4.1)) 
  p <- plot(mat, xlab = "Month", ylab = "Year", main = region_name_, key = NULL,
       col=inferno(10), breaks=c(0, 0.1, 0.2, 0.3, 0.4, 
                                 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  
}

region_names <- c("Cape Cod Bay",
                  "Great South Channel",
                  "Massachusetts Bay",
                  "Jordan Basin",
                  "Wilkinson Basin")

regions <- c("CCB", "GSC", "MB", "JB", "WB")

thresholds <- c(10000)

md$thresh_10 <- if_else(md$abund > thresholds[1], 1, 0)

par(mfrow = c(1, 5), mar=c(3, 4, 1.5, 1.5), cex.axis = 1.5, las = 1) 
for (j in 1:length(thresholds)) {
  for (i in 1:length(regions)) {
    mat <- matrix(nrow = length(years), ncol = 12)
    region_md <- md %>% dplyr::filter(region == regions[i])
    
    for (k in 1:length(years)) {
      for (l in 1:12) {
        year_ <- years[k]
        month_md <- region_md %>% filter(year == year_, month == l)
        month_md$abund <- if_else(month_md$abund > thresholds[j], 1, 0)
      
        mat[k,l] <- mean(month_md$abund)
      }
    }
    
    mat[is.nan(mat)] <- NA
    rownames(mat) <- years
  
    plot(mat, xlab = "", ylab = "", main = region_names[i], key = NULL,
              col=inferno(11), breaks=c(0, 0.1, 0.2, 0.3, 0.4,
                                        0.5, 0.6, 0.7, 0.8,
                                        0.9, 1.0)
         )
    
    corr_dat <- as.data.frame(pivot_longer(as.data.frame(mat), cols = 1:12))
    corr_dat[is.na(corr_dat)] <- 0
    corr_dat$month <- rep(1:12, 18)
    corr_dat$year <- rep(2000:2017, each = 12, length.out = 216)
    
    corr_dat <- corr_dat %>%
      dplyr::group_by(month) %>%
      dplyr::summarise(cor = stats::cor.test(value, year, na.action = na.omit())$estimate,
                       sig = stats::cor.test(value, year, na.action = na.omit())$p.value)
    
    print(i)
    print(corr_dat)
    
   }
}


par(mfrow = c(1, 5), mar=c(3, 4, 1.5, 1.5), cex.axis = 1.5, las = 1) 
for (j in 1:length(thresholds)) {
  for (i in 1:length(regions)) {
    mat <- matrix(nrow = length(years), ncol = 12)
    
    for (k in 1:length(years)) {
      for (l in 1:12) {
        print(paste(years[k], ":", l))
        rf_proj_raster <- raster(file.path(fp_out, species, version, "Biomod", "Projections", paste0(species, version), paste0('proj_RF_',  years[k], "_", l), paste0('proj_RF_',  years[k], "_", l, "_", species, version, '.grd'))) %>%
          `/`(1000)
        
        crs(rf_proj_raster) <- '+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62'
        
        # Save projection raster as data frame with xy coords and no NAs
        rf_proj_df <- as.data.frame(rf_proj_raster, xy = TRUE, na.rm = TRUE)
        
        # Assign column names
        names(rf_proj_df) <- c('x', 'y', "pred")
        
        rf_proj_df$x <- round(rf_proj_df$x, 3)
        rf_proj_df$y <- round(rf_proj_df$y, 3)
        
        proj_sf <- st_as_sf(rf_proj_df, coords = c("x", "y"))
        
        iRB <- st_intersects(RB, proj_sf, sparse = FALSE)
        iCCB <- st_intersects(CCB, proj_sf, sparse = FALSE)
        iBOF <- st_intersects(BOF, proj_sf, sparse = FALSE)
        iWB <- st_intersects(WB, proj_sf, sparse = FALSE)
        iJB <- st_intersects(JB, proj_sf, sparse = FALSE)
        iGSC <- st_intersects(GSC, proj_sf, sparse = FALSE)
        iMB <- st_intersects(MB, proj_sf, sparse = FALSE)
        
        poly_names <- rep("", nrow(rf_proj_df))
        poly_names[iMB] <- "MB"
        poly_names[iCCB] <- "CCB"
        poly_names[iWB] <- "WB"
        poly_names[iJB] <- "JB"
        poly_names[iGSC] <- "GSC"
        
        rf_proj_df <- rf_proj_df %>%
          dplyr::mutate(region = poly_names)
        
        rf_proj_df$region[rf_proj_df$region == ""] <- NA
        
        
        region_proj <- rf_proj_df %>% dplyr::filter(region == regions[i])
        
        mat[k,l] <- mean(region_proj$pred)
      }
    }
    
    mat[is.nan(mat)] <- NA
    rownames(mat) <- years
    
    plot(mat, xlab = "", ylab = "", main = region_names[i], key = NULL,
         col=inferno(11), breaks=c(0, 0.1, 0.2, 0.3, 0.4,
                                   0.5, 0.6, 0.7, 0.8,
                                   0.9, 1.0)
    )
    
    corr_dat <- as.data.frame(pivot_longer(as.data.frame(mat), cols = 1:12))
    corr_dat[is.na(corr_dat)] <- 0
    corr_dat$month <- rep(1:12, 18)
    corr_dat$year <- rep(2000:2017, each = 12, length.out = 216)
    
    corr_dat <- corr_dat %>%
      dplyr::group_by(month) %>%
      dplyr::summarise(cor = stats::cor.test(value, year, na.action = na.omit())$estimate,
                       sig = stats::cor.test(value, year, na.action = na.omit())$p.value)
    
    print(i)
    print(corr_dat)
    
    # regress_mat <- as.data.frame(pivot_longer(as.data.frame(mat), cols = 1:12))
    # regress_mat[is.na(regress_mat)] <- 0
    # regress_mat$month <- rep(1:12, 18)
    # regress_mat$year <- rep(2000:2017, each = 12, length.out = 216)
    # 
    # 
    # regress_mat <- regress_mat %>% 
    #   group_by(year) %>%
    #   mutate(max = which.max(value)) %>%
    #   summarise(max = mean(max))
    # 
    # lin_regress <- lm(max ~ year, data = regress_mat)
    # summary(lin_regress)
    
    
    
  }
}

version = "v0.4.1"
fp_proj = file.path(fp_out, species, version, "Biomod", "Projections")

proj_10 <- raster::raster(file.path(fp_proj, paste0(species, version), "proj_RF_1", paste0("proj_RF_1_", species, version, ".grd")))/1000
names(proj_10) <- c("1")
for (i in 2:12) {
  proj_temp <- raster::raster(file.path(fp_proj, paste0(species, version), paste0("proj_RF_", i), paste0("proj_RF_", i, "_", species, version, ".grd")))/1000
  names(proj_temp) <- c(i)
  proj_10 <- raster::stack(proj_10, proj_temp)
}

proj_10 <- proj_10 %>% as.data.frame(xy = TRUE)
names(proj_10) <- c("lon", "lat", 1:12)

proj_sf <- st_as_sf(proj_10, coords = c("lon", "lat"))

iRB <- st_intersects(RB, proj_sf, sparse = FALSE)
iCCB <- st_intersects(CCB, proj_sf, sparse = FALSE)
iBOF <- st_intersects(BOF, proj_sf, sparse = FALSE)
iWB <- st_intersects(WB, proj_sf, sparse = FALSE)
iJB <- st_intersects(JB, proj_sf, sparse = FALSE)
iGSC <- st_intersects(GSC, proj_sf, sparse = FALSE)
iMB <- st_intersects(MB, proj_sf, sparse = FALSE)

poly_names <- rep("", nrow(proj_10))
poly_names[iMB] <- "MB"
poly_names[iCCB] <- "CCB"
poly_names[iWB] <- "WB"
poly_names[iJB] <- "JB"
poly_names[iGSC] <- "GSC"

proj_10 <- proj_10 %>%
  dplyr::mutate(region = poly_names)

proj_10$region[proj_10$region == ""] <- NA

version = "v0.6.5"
fp_proj = file.path(fp_out, species, version, "Biomod", "Projections")

proj_40 <- raster::raster(file.path(fp_proj, paste0(species, version), "proj_RF_1", paste0("proj_RF_1_", species, version, ".grd")))/1000
names(proj_40) <- c("1")
for (i in 2:12) {
  proj_temp <- raster::raster(file.path(fp_proj, paste0(species, version), paste0("proj_RF_", i), paste0("proj_RF_", i, "_", species, version, ".grd")))/1000
  names(proj_temp) <- c(i)
  proj_40 <- stack(proj_40, proj_temp)
}

proj_40 <- proj_40 %>% as.data.frame(xy = TRUE)
names(proj_40) <- c("lon", "lat", 1:12)

proj_sf <- st_as_sf(proj_40, coords = c("lon", "lat"))

iRB <- st_intersects(RB, proj_sf, sparse = FALSE)
iCCB <- st_intersects(CCB, proj_sf, sparse = FALSE)
iBOF <- st_intersects(BOF, proj_sf, sparse = FALSE)
iWB <- st_intersects(WB, proj_sf, sparse = FALSE)
iJB <- st_intersects(JB, proj_sf, sparse = FALSE)
iGSC <- st_intersects(GSC, proj_sf, sparse = FALSE)
iMB <- st_intersects(MB, proj_sf, sparse = FALSE)

poly_names <- rep("", nrow(proj_40))
poly_names[iMB] <- "MB"
poly_names[iCCB] <- "CCB"
poly_names[iWB] <- "WB"
poly_names[iJB] <- "JB"
poly_names[iGSC] <- "GSC"

proj_40 <- proj_40 %>%
  dplyr::mutate(region = poly_names)

proj_40$region[proj_40$region == ""] <- NA


regions <- c("CCB", "GSC", "MB", "JB", "WB")

thresholds <- c(1000, 4000, 10000, 40000)

par(mfrow = c(2, 5), mar=c(3, 4, 1.5, 1.5), cex.axis = 1.5, las = 1) 
for (j in 1:length(thresholds)) {
  for (i in 1:length(regions)) {
    mat <- matrix(nrow = 1, ncol = 12)
    region_proj <- proj %>% dplyr::filter(region == regions[i])
    
    for (l in 1:12) {
      month_proj <- region_proj %>% dplyr::select(`l`)
      mat[1,l] <- mean(month_proj)
    }
  }
    
  mat[is.nan(mat)] <- NA
  
  plot(mat, xlab = "", ylab = "", xlim = 1:12, main = region_names[i], key = NULL,
       col=inferno(10), breaks=c(0, 0.1, 0.2, 0.3, 0.4, 
                                   0.5, 0.6, 0.7, 0.8, 0.9, 1)
  )
}

md$thresh_1 <- if_else(md$abund > thresholds[1], 1, 0)
md$thresh_4 <- if_else(md$abund > thresholds[2], 1, 0)
md$thresh_10 <- if_else(md$abund > thresholds[3], 1, 0)
md$thresh_40 <- if_else(md$abund > thresholds[4], 1, 0)

proj_long_10 <- proj_10 %>% pivot_longer(cols = 3:14, names_to = c("month"), values_to = c("proj"))
proj_long_10$lon <- round(proj_long_10$lon, 2)
proj_long_10$lat <- round(proj_long_10$lat, 2)

proj_long_40 <- proj_40 %>% pivot_longer(cols = 3:14, names_to = c("month"), values_to = c("proj"))
proj_long_40$lon <- round(proj_long_40$lon, 2)
proj_long_40$lat <- round(proj_long_40$lat, 2)

proj_long_10$month <- as.double(proj_long_10$month)
proj_long_40$month <- as.double(proj_long_40$month)

joint_10 <- dplyr::left_join(x = md, y = proj_long_10, by = c("lat", "lon", "region", "month"))
  
joint_10 <- joint_10 %>% 
         dplyr::group_by(month) %>%
         dplyr::summarise(mean_thresh_10 = mean(thresh_10, na.rm = TRUE),
                          mean_proj = mean(proj, na.rm = TRUE, nan.remove = TRUE))

joint_40 <- dplyr::left_join(x = md, y = proj_long_40, by = c("lat", "lon", "region", "month"))

joint_40 <- joint_40 %>% 
  dplyr::group_by(month, year) %>%
  dplyr::summarise(mean_thresh_40 = mean(thresh_1, na.rm = TRUE, nan.remove = TRUE),
                   mean_proj = mean(proj, na.rm = TRUE, nan.remove = TRUE))

linear_model <- lm(mean_thresh_40 ~ mean_proj, data = joint_40)
summary(linear_model)

ggplot(joint_40, aes(x = mean_thresh_40, y = mean_proj)) +
  geom_point(cex = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  labs(x = "Actual", y = "Projection") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 20))



