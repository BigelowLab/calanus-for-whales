library(plyr)
library(dplyr)
library(ggplot2)
library(viridis)


months <- c("January", "February", "March",
            "April", "May", "June",
            "July", "August", "September",
            "October", "November", "December")

# -------- Load world map data --------
worldmap <- ggplot2::map_data("world")

biomod_dataset = "ECOMON"

# Create Grieve et al. 2017 style plots
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

thresholds <- c(1000, 4000, 10000, 40000)

md$thresh_1 <- if_else(md$abund > thresholds[1], 1, 0)
md$thresh_4 <- if_else(md$abund > thresholds[2], 1, 0)
md$thresh_10 <- if_else(md$abund > thresholds[3], 1, 0)
md$thresh_40 <- if_else(md$abund > thresholds[4], 1, 0)

# 1,000 full data
for (i in 1:12) {
  
  ggplot(md %>% filter(month == i), aes(x = lon, y = lat, col = thresh_1)) +
    geom_point(cex = 1) +
    scale_color_viridis(option = "inferno", limits = c(0,1)) +
    labs(x = "",
         y = "",
         col = "Proportion") +
    ggtitle(paste0(months[i], " 1,000 Full")) +
    # Add world map data
    geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
    coord_quickmap(xlim = c(-76, -64),
                   ylim = c(35, 47),
                   expand = TRUE) +
    # Remove grid lines
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank(),
          plot.title = element_text(size=52)) +
    ggsave(file.path(DIR, "Figures", "Grieve_et_al_2017_plots", paste0("grieve_full_1000_", i, ".png")))
  
}

# 4,000 full data
for (i in 1:12) {
  
  ggplot(md %>% filter(month == i), aes(x = lon, y = lat, col = thresh_4)) +
    geom_point(cex = 1) +
    scale_color_viridis(option = "inferno", limits = c(0,1)) +
    labs(x = "",
         y = "",
         col = "Proportion") +
    ggtitle(paste0(months[i], " 4,000 Full")) +
    # Add world map data
    geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
    coord_quickmap(xlim = c(-76, -64),
                   ylim = c(35, 47),
                   expand = TRUE) +
    # Remove grid lines
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank(),
          plot.title = element_text(size=52)) +
    ggsave(file.path(DIR, "Figures", "Grieve_et_al_2017_plots", paste0("grieve_full_4000_", i, ".png")))
  
}

# 10,000 full data
for (i in 1:12) {
  
  ggplot(md %>% filter(month == i), aes(x = lon, y = lat, col = as.factor(thresh_10))) +
    geom_point(cex = 2, alpha = 0.7) +
    scale_color_viridis_d(option = "viridis") +
    labs(x = "",
         y = "") +
    ggtitle(months[i]) +
    # Add world map data
    geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
    coord_quickmap(xlim = c(-76, -64),
                   ylim = c(35, 47),
                   expand = TRUE) +
    # Remove grid lines
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank(),
          plot.title = element_text(size=52),
          legend.position = "none") +
    ggsave(file.path(DIR, "Figures", "Grieve_et_al_2017_plots", paste0(species, "_grieve_full_10000_", i, ".png")))
  
}

# 40,000 full data
for (i in 1:12) {
  
  ggplot(md %>% filter(month == i), aes(x = lon, y = lat, col = thresh_40)) +
    geom_point(cex = 1) +
    scale_color_viridis(option = "inferno", limits = c(0,1)) +
    labs(x = "",
         y = "",
         col = "Proportion") +
    ggtitle(paste0(months[i], " 40,000 Full")) +
    # Add world map data
    geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
    coord_quickmap(xlim = c(-76, -64),
                   ylim = c(35, 47),
                   expand = TRUE) +
    # Remove grid lines
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank(),
          plot.title = element_text(size=52)) +
    ggsave(file.path(DIR, "Figures", "Grieve_et_al_2017_plots", paste0("grieve_full_40000_", i, ".png")))
  
}

biomod_dataset = "ECOMON_STAGED"

# Create Grieve et al. 2017 style plots
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

thresholds <- c(1000, 4000, 10000, 40000)

md$thresh_1 <- if_else(md$abund > thresholds[1], 1, 0)
md$thresh_4 <- if_else(md$abund > thresholds[2], 1, 0)
md$thresh_10 <- if_else(md$abund > thresholds[3], 1, 0)
md$thresh_40 <- if_else(md$abund > thresholds[4], 1, 0)

# 1,000 full data
grieve <- md %>% 
  dplyr::group_by(month, lat, lon) %>%
  dplyr::summarise(mean_thresh = mean(thresh_1, na.rm = TRUE, nan.remove = TRUE))

for (i in 1:12) {
  
  ggplot(grieve %>% filter(month == i), aes(x = lon, y = lat, col = mean_thresh)) +
    geom_point(cex = 1) +
    scale_color_viridis(option = "inferno", limits = c(0,1)) +
    labs(x = "",
         y = "",
         col = "Proportion") +
    ggtitle(paste0(months[i], " 1,000 Staged")) +
    # Add world map data
    geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
    coord_quickmap(xlim = c(-76, -64),
                   ylim = c(35, 47),
                   expand = TRUE) +
    # Remove grid lines
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank(),
          plot.title = element_text(size=52)) +
    ggsave(file.path(DIR, "Figures", "Grieve_et_al_2017_plots", paste0("grieve_staged_1000_", i, ".png")))
  
}

# 4,000 full data
grieve <- md %>% 
  dplyr::group_by(month, lat, lon) %>%
  dplyr::summarise(mean_thresh = mean(thresh_4, na.rm = TRUE, nan.remove = TRUE))

for (i in 1:12) {
  
  ggplot(grieve %>% filter(month == i), aes(x = lon, y = lat, col = mean_thresh)) +
    geom_point(cex = 1) +
    scale_color_viridis(option = "inferno", limits = c(0,1)) +
    labs(x = "",
         y = "",
         col = "Proportion") +
    ggtitle(paste0(months[i], " 4,000 Staged")) +
    # Add world map data
    geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
    coord_quickmap(xlim = c(-76, -64),
                   ylim = c(35, 47),
                   expand = TRUE) +
    # Remove grid lines
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank(),
          plot.title = element_text(size=52)) +
    ggsave(file.path(DIR, "Figures", "Grieve_et_al_2017_plots", paste0("grieve_staged_4000_", i, ".png")))
  
}

# 10,000 full data
grieve <- md %>% 
  dplyr::group_by(month, lat, lon) %>%
  dplyr::summarise(mean_thresh = mean(thresh_10, na.rm = TRUE, nan.remove = TRUE))

for (i in 1:12) {
  
  ggplot(grieve %>% filter(month == i), aes(x = lon, y = lat, col = mean_thresh)) +
    geom_point(cex = 1) +
    scale_color_viridis(option = "inferno", limits = c(0,1)) +
    labs(x = "",
         y = "",
         col = "Proportion") +
    ggtitle(paste0(months[i], " 10,000 Staged")) +
    # Add world map data
    geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
    coord_quickmap(xlim = c(-76, -64),
                   ylim = c(35, 47),
                   expand = TRUE) +
    # Remove grid lines
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank(),
          plot.title = element_text(size=52)) +
    ggsave(file.path(DIR, "Figures", "Grieve_et_al_2017_plots", paste0("grieve_staged_10000_", i, ".png")))
  
}

# 40,000 full data
grieve <- md %>% 
  dplyr::group_by(month, lat, lon) %>%
  dplyr::summarise(mean_thresh = mean(thresh_4, na.rm = TRUE, nan.remove = TRUE))

for (i in 1:12) {
  
  ggplot(grieve %>% filter(month == i), aes(x = lon, y = lat, col = mean_thresh)) +
    geom_point(cex = 1) +
    scale_color_viridis(option = "inferno", limits = c(0,1)) +
    labs(x = "",
         y = "",
         col = "Proportion") +
    ggtitle(paste0(months[i], " 40,000 Staged")) +
    # Add world map data
    geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
    coord_quickmap(xlim = c(-76, -64),
                   ylim = c(35, 47),
                   expand = TRUE) +
    # Remove grid lines
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank(),
          plot.title = element_text(size=52)) +
    ggsave(file.path(DIR, "Figures", "Grieve_et_al_2017_plots", paste0("grieve_staged_40000_", i, ".png")))
  
}



