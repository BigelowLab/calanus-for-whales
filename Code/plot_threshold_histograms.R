
library(dplyr)
library(readr)

# -------- Initialize filepaths --------
# Data filepath
fp_data <- "./calanus_data/Data"
# Databases filepath
fp_db <- file.path(fp_data, "Databases")
  
# ---- Load color-blind friendly palette ----
# Source: https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/#use-a-colorblind-friendly-palette
cbp <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
# --- Load data ----
# md <- load_md(version, fp_md, biomod_dataset, fp_covars, env_covars, 
#               years, fp_out, species)
#               
select_dataset <- "ECOMON"

if (select_dataset == "ECOMON") {
  md <- read_csv(file.path(fp_db, "zooplankton_database.csv")) %>%
    dplyr::filter(dataset == select_dataset, year %in% 2000:2017) %>%
    #dplyr::mutate(abund = cfin_CV_VI)
    #dplyr::mutate(abund = ctyp_total)
    dplyr::mutate(abund = pseudo_total)
} else if (select_dataset == "ECOMON_STAGED") {
  md <- read_csv(file.path(fp_db, "zooplankton_database.csv")) %>%
    dplyr::filter(dataset == select_dataset, year %in% 2000:2017) %>%
    dplyr::mutate(abund = cfin_CIV + cfin_CV + cfin_CVI)
}
  
# Define lower limit (thickness of calanus layer)
lower <- 1
# Define upper limit (thickness of calanus layer)
upper <- 20

# Cape Cod Bay
# Fortune et al. 2013
fortune_ccb_lower <- 14778 * lower
fortune_ccb_upper <- 14778 * upper
# Mayo and Marx 1990
mayo_lower <- 1000 * lower
mayo_upper <- 1000 * upper

# Great South Channel
# Beardley et al. 1996
beardsley_lower <- 8700 * lower
beardsley_upper <- 41000 * upper
# Kenney et al. 1986
kenney_lower <- 300000 * lower
kenney_upper <- 1000000 * upper
# Wishner et al. 1988
wishner88_lower <- 41600 * lower
wishner88_upper <- 41600 * upper
# Wishner et al. 1995
wishner95_lower <- 9749 * lower
wishner95_upper <- 9749 * upper

# Bay of Fundy
# Baumgartner and Mate 2003
baumgartner_mate_lower <- 3600 * lower
baumgartner_mate_upper <- 3600 * upper
# Fortune et al. 2013
fortune_bof_lower <- 6618 * lower
fortune_bof_upper <- 6618 * upper
# Michaud and Taggart 2007
michaud_lower <- 900 * lower
michaud_upper <- 900 * upper
# Murison and Gaskin 1988
murison_lower <- 832 * lower
murison_upper <- 1070 * upper
# Woodley and Gaskin 1996
woodley_lower <- 1139 * lower
woodley_upper <- 1139 * upper

# Gulf of Maine
# Baumgartner et al. 2017
baumgartner_lower <- 14900 * lower
baumgartner_upper <- 14900 * upper
# Record et al. 2019
record_lower <- 35000 
record_upper <- 45000

base <- 2500

hist(log10(md$abund), main = "EcoMon Logged Abundances", 
     xlab = expression(paste("Logged Abundance (m"^"-2", ")")), xlim = c(0,6), ylim = c(0,5000),
     col = "white", border = "black")
# ---- Add vertical line for 1,000 ----
abline(v = log10(1000), col = "#C0C0C0", lwd = 5, lty = 1)
# ---- Add vertical line for 4,000 ----
abline(v = log10(4000), col = "#C0C0C0", lwd = 5, lty = 1)
# ---- Add vertical line for 10,000 ----
abline(v = log10(10000), col = "#C0C0C0", lwd = 5, lty = 1)
# ---- Add vertical line for 40,000 ----
abline(v = log10(40000), col = "#C0C0C0", lwd = 5, lty = 1)
# ---- Add horizontal line for Cape Cod Bay ----
lines(x = c(log10(mayo_lower), log10(mayo_upper)), 
      y = c(base, base), col = cbp[8], lwd = 5, lty = 1)
lines(x = c(log10(fortune_ccb_lower), log10(fortune_ccb_upper)), 
      y = c(base + 250, base + 250), col = cbp[8], lwd = 5, lty = 1)
# ---- Add horizontal line for Bay of Fundy ----
lines(x = c(log10(murison_lower), log10(murison_upper)), 
      y = c(base + 500, base + 500), col = cbp[3], lwd = 5, lty = 1)
lines(x = c(log10(michaud_lower), log10(michaud_upper)), 
      y = c(base + 750, base + 750), col = cbp[3], lwd = 5, lty = 1)
lines(x = c(log10(woodley_lower), log10(woodley_upper)), 
      y = c(base + 1000, base + 1000), col = cbp[3], lwd = 5, lty = 1)
lines(x = c(log10(baumgartner_mate_lower), log10(baumgartner_mate_upper)), 
      y = c(base + 1250, base + 1250), col = cbp[3], lwd = 5, lty = 1)
lines(x = c(log10(fortune_bof_lower), log10(fortune_bof_upper)), 
      y = c(base + 1500, base + 1500), col = cbp[3], lwd = 5, lty = 1)
# ---- Add horizontal line for Gulf of Maine ----
lines(x = c(log10(baumgartner_lower), log10(baumgartner_upper)), 
      y = c(base + 1750, base + 1750), col = cbp[5], lwd = 5, lty = 1)
lines(x = c(log10(record_lower), log10(record_upper)), 
      y = c(base + 2000, base + 2000), col = cbp[5], lwd = 5, lty = 1)
# ---- Add horizontal line for Great South Channel ----
lines(x = c(log10(beardsley_lower), log10(beardsley_upper)), 
      y = c(base + 2250, base + 2250), col = cbp[2], lwd = 5, lty = 1)
lines(x = c(log10(wishner95_lower), log10(wishner95_upper)), 
      y = c(base + 2500, base + 2500), col = cbp[2], lwd = 5, lty = 1)
lines(x = c(log10(wishner88_lower), log10(wishner88_upper)), 
      y = c(base + 2750, base + 2750), col = cbp[2], lwd = 5, lty = 1)
lines(x = c(log10(kenney_lower), log10(kenney_upper)), 
      y = c(base + 3000, base + 3000), col = cbp[2], lwd = 5, lty = 1)

legend(6, 5000, c("Great South Channel",
                  "Gulf of Maine", 
                  "Bay of Fundy",
                  "Cape Cod Bay"),
       title = "Region", col=c(cbp[2], cbp[5],
                                cbp[3], cbp[8]), 
       lwd=5, lty = 1, cex = 0.7)


