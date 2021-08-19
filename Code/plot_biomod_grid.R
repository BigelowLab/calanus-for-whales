library(png)
library(gridExtra)
library(grid)

version = "v0.7.5"

fp <- file.path(fp_out, species, version, "Biomod", "Plots")

p1 <- readPNG(file.path(fp, "rf_proj_1.png"))
p2 <- readPNG(file.path(fp, "rf_proj_2.png"))
p3 <- readPNG(file.path(fp, "rf_proj_3.png"))
p4 <- readPNG(file.path(fp, "rf_proj_4.png"))
p5 <- readPNG(file.path(fp, "rf_proj_5.png"))
p6 <- readPNG(file.path(fp, "rf_proj_6.png"))
p7 <- readPNG(file.path(fp, "rf_proj_7.png"))
p8 <- readPNG(file.path(fp, "rf_proj_8.png"))
p9 <- readPNG(file.path(fp, "rf_proj_9.png"))
p10 <- readPNG(file.path(fp, "rf_proj_10.png"))
p11 <- readPNG(file.path(fp, "rf_proj_11.png"))
p12 <- readPNG(file.path(fp, "rf_proj_12.png"))


par(mar=c(0.5, 0.5, 0.5, 0.5)) 
grid.arrange(rasterGrob(p1), rasterGrob(p2), rasterGrob(p3),
             rasterGrob(p4), rasterGrob(p5), rasterGrob(p6),
             rasterGrob(p7), rasterGrob(p8), rasterGrob(p9),
             rasterGrob(p10), rasterGrob(p11), rasterGrob(p12),
             ncol = 4)

version = "v0.4.4"
fp <- file.path(fp_out, species, version, "Biomod", "Plots")

for (i in 1:12) {
  i=12
  p1 <- readPNG(file.path(fp, paste0("rf_proj_2000_", i, ".png")))
  p2 <- readPNG(file.path(fp, paste0("rf_proj_2001_", i, ".png")))
  p3 <- readPNG(file.path(fp, paste0("rf_proj_2002_", i, ".png")))
  p4 <- readPNG(file.path(fp, paste0("rf_proj_2003_", i, ".png")))
  p5 <- readPNG(file.path(fp, paste0("rf_proj_2004_", i, ".png")))
  p6 <- readPNG(file.path(fp, paste0("rf_proj_2005_", i, ".png")))
  p7 <- readPNG(file.path(fp, paste0("rf_proj_2006_", i, ".png")))
  p8 <- readPNG(file.path(fp, paste0("rf_proj_2007_", i, ".png")))
  p9 <- readPNG(file.path(fp, paste0("rf_proj_2008_", i, ".png")))
  p10 <- readPNG(file.path(fp, paste0("rf_proj_2009_", i, ".png")))
  p11 <- readPNG(file.path(fp, paste0("rf_proj_2010_", i, ".png")))
  p12 <- readPNG(file.path(fp, paste0("rf_proj_2011_", i, ".png")))
  p13 <- readPNG(file.path(fp, paste0("rf_proj_2012_", i, ".png")))
  p14 <- readPNG(file.path(fp, paste0("rf_proj_2013_", i, ".png")))
  p15 <- readPNG(file.path(fp, paste0("rf_proj_2014_", i, ".png")))
  p16 <- readPNG(file.path(fp, paste0("rf_proj_2015_", i, ".png")))
  p17 <- readPNG(file.path(fp, paste0("rf_proj_2016_", i, ".png")))
  p18 <- readPNG(file.path(fp, paste0("rf_proj_2017_", i, ".png")))
  
  par(mar=c(0.1, 0.1, 0.1, 0.1)) 
  grid.arrange(rasterGrob(p1), rasterGrob(p2), rasterGrob(p3),
               rasterGrob(p4), rasterGrob(p5), rasterGrob(p6),
               rasterGrob(p7), rasterGrob(p8), rasterGrob(p9),
               rasterGrob(p10), rasterGrob(p11), rasterGrob(p12),
               rasterGrob(p13), rasterGrob(p14), rasterGrob(p15),
               rasterGrob(p16), rasterGrob(p17), rasterGrob(p18),
               ncol = 6)
  
}

# Grieve et al plots
fp <- file.path(DIR, "Figures", "Grieve_et_al_2017_plots")

p1 <- readPNG(file.path(fp, paste0(species, "_grieve_full_10000_1.png")))
p2 <- readPNG(file.path(fp, paste0(species, "_grieve_full_10000_2.png")))
p3 <- readPNG(file.path(fp, paste0(species, "_grieve_full_10000_3.png")))
p4 <- readPNG(file.path(fp, paste0(species, "_grieve_full_10000_4.png")))
p5 <- readPNG(file.path(fp, paste0(species, "_grieve_full_10000_5.png")))
p6 <- readPNG(file.path(fp, paste0(species, "_grieve_full_10000_6.png")))
p7 <- readPNG(file.path(fp, paste0(species, "_grieve_full_10000_7.png")))
p8 <- readPNG(file.path(fp, paste0(species, "_grieve_full_10000_8.png")))
p9 <- readPNG(file.path(fp, paste0(species, "_grieve_full_10000_9.png")))
p10 <- readPNG(file.path(fp, paste0(species, "_grieve_full_10000_10.png")))
p11 <- readPNG(file.path(fp, paste0(species, "_grieve_full_10000_11.png")))
p12 <- readPNG(file.path(fp, paste0(species, "_grieve_full_10000_12.png")))


par(mar=c(0.1, 0.1, 0.1, 0.1)) 
grid.arrange(rasterGrob(p1), rasterGrob(p2), rasterGrob(p3),
             rasterGrob(p4), rasterGrob(p5), rasterGrob(p6),
             rasterGrob(p7), rasterGrob(p8), rasterGrob(p9),
             rasterGrob(p10), rasterGrob(p11), rasterGrob(p12),
             ncol = 4)


