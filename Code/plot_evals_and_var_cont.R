library(readr)
library(tidyr)
library(dplyr)

version = "v0.7.2"
fp_evals = file.path(fp_out, species, version, "Biomod", "Evals")

#Color blind friendly pallette. Source: https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/#use-a-colorblind-friendly-palette
cbp <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


full_df_init <- read_csv(file.path(fp_evals, 'evals.csv')) %>%
  dplyr::filter(Eval.metric == "TSS", grepl('RUN', Model.name )) %>%
  dplyr::mutate(month = 1, model_id = rep(seq(1:3),10), EM = FALSE)

full_df_init <- full_df_init %>%
  dplyr::group_by(model_id) %>%
  dplyr::mutate(min = min(Testing.data), max = max(Testing.data),
                se = sd(Testing.data)/sqrt(length(full_df_init)))

#full_df <- full_df_init

full_df <- aggregate(full_df_init[, 4], list(full_df_init$model_id), mean) %>%
  dplyr::mutate(month = 1, model_id = seq(1:3), Model.name = c("GAM", "GBM", "RF"), EM = FALSE,
               min = unique(full_df_init$min), max = unique(full_df_init$max),
               se = unique(full_df_init$se)) %>%
  dplyr::select(Testing.data, month, model_id, Model.name, EM, min, max, se)

EM <- read_csv(file.path(fp_evals, 'ensemble_evals_1.csv')) %>%
  dplyr::filter(Eval.metric == "ROC", grepl('EMmean', Model.name )) %>%
  dplyr::mutate(month = 1, model_id = 4, EM = TRUE, min = Testing.data, max = Testing.data) %>%
  dplyr::select(Testing.data, month, model_id, Model.name, EM, min, max)

full_df <- rbind(full_df, EM)

months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

for (i in c(2:6, 8:11)) {
  df_init <- read_csv(file.path(fp_evals, paste0('evals_', i, '.csv'))) %>%
    filter(Eval.metric == "ROC", grepl('RUN', Model.name )) %>%
    mutate(month = i, model_id = rep(seq(1:3),10), EM = FALSE) %>%
    dplyr::group_by(model_id) %>%
    mutate(min = min(Testing.data), max = max(Testing.data)) 
  
  #df <- df_init
  print(i)
  if (i != 7 & i != 12) {
    df <- aggregate(df_init[, 4], list(df_init$model_id), mean) %>%
      mutate(month = i, model_id = seq(1:3), Model.name = c("GAM", "GBM", "RF"), EM = FALSE,
             min = unique(df_init$min), max = unique(df_init$max)) %>%
      dplyr::select(Testing.data, month, model_id, Model.name, EM, min, max)
  }
  
  full_df <- rbind(full_df, df)
  
  EM <- read_csv(file.path(fp_evals, paste('ensemble_evals_', i, '.csv', sep = ""))) %>%
    filter(Eval.metric == "ROC", grepl('EMmean', Model.name )) %>%
    mutate(month = i, model_id = 4, EM = TRUE, min = Testing.data, max = Testing.data) %>%
    dplyr::select(Testing.data, month, model_id, Model.name, EM, min, max)
  
  full_df <- rbind(full_df, EM)
  
}

full_df$Model.name[full_df$model_id == 1] <- "BRT"
full_df$Model.name[full_df$model_id == 2] <- "GAM"
full_df$Model.name[full_df$model_id == 3] <- "RF"
full_df$Model.name[full_df$model_id == 4] <- "Ensemble"

# Plot ensemble AUC
plot(x = full_df$month[full_df$Model.name == "RF"], 
     y = full_df$Testing.data[full_df$Model.name == "RF"], type = "l", col = cbp[3], lty = 1,
     xlim = c(1,12), ylim = c(0.6,1.0), ylab = "AUC", xlab = "Month", axes = FALSE)
# Add x axis
axis(1, at = c(1,2,3,4,5,6,7,8,9,10,11,12), labels = months)
# Add y axis
axis(2, at = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))
# Plot mean ANN AUC
lines(x = full_df$month[full_df$Model.name == "RF"], 
      y = full_df$Testing.data[full_df$Model.name == "RF"], type = "l", col = cbp[2])
# Plot ANN range
polygon(c(full_df$month[full_df$Model.name == "RF"],rev(full_df$month[full_df$Model.name == "RF"])),
        c(full_df$min[full_df$Model.name == "RF"],rev(full_df$max[full_df$Model.name == "RF"])),
        col = alpha(cbp[3], 0.2), border = FALSE)
# Plot mean BRT AUC
lines(x = full_df$month[full_df$Model.name == "BRT"], 
      y = full_df$Testing.data[full_df$Model.name == "BRT"], type = "l", col = cbp[8])
# Plot BRT range
polygon(c(full_df$month[full_df$Model.name == "BRT"],rev(full_df$month[full_df$Model.name == "BRT"])),
        c(full_df$min[full_df$Model.name == "BRT"],rev(full_df$max[full_df$Model.name == "BRT"])),
        col = alpha(cbp[8], 0.2), border = FALSE)
# Plot mean GAM AUC
lines(x = full_df$month[full_df$Model.name == "GAM"], 
      y = full_df$Testing.data[full_df$Model.name == "GAM"], type = "l", col = cbp[6])
# Plot GAM range
polygon(c(full_df$month[full_df$Model.name == "GAM"],rev(full_df$month[full_df$Model.name == "GAM"])),
        c(full_df$min[full_df$Model.name == "GAM"],rev(full_df$max[full_df$Model.name == "GAM"])),
        col = alpha(cbp[6], 0.2), border = FALSE)
# Add legend
legend(x = 6.8, y = 0.45, legend = c("Ensemble", "RF", "BRT", "GAM"),
       col = cbp[c(3,2,8,6)], lty = c(2, 1, 1, 1), cex = 0.8, box.col = "white")

version = "v0.5.9"
fp_evals = file.path(fp_out, species, version, "Biomod", "Evals")

env_covars <- c("wind", "int_chl", "sst", "sst_grad", "jday",
                "uv_grad", "bat", "slope", "bots", "bott")

df_var <- read_csv(file.path(fp_evals, 'var_importance_1.csv')) %>%
  dplyr::select(contains('RUN'))

new_df_rf <- df_var %>% dplyr::select(contains("RF"))
new_df_rf <- new_df_rf %>% gather() %>% mutate(Var = rep(env_covars,10),
                                                 Month = 1)
names(new_df_rf) <- c("Model", "Perc.Cont", "Var", "Month") 
new_df_rf <- new_df_rf %>% mutate(Model = "RF", 
                                    Perc.Cont.Temp = rep(rowMeans(df_var %>% 
                                                                    dplyr::select(
                                                                      contains("RF")), 
                                                                  na.rm = TRUE), 10)) %>%
  mutate(Perc.Cont.Norm = Perc.Cont/sum(unique(Perc.Cont.Temp))/10)

sum(new_df_rf$Perc.Cont.Norm)

#new_df_rf <- aggregate(new_df_rf[, 6], list(new_df_rf$Var), mean)

library(ggplot2)
ggplot(new_df_rf, aes(x = 1, y = Perc.Cont.Norm, fill = Group.1)) +
  geom_col() +
  scale_fill_viridis_d() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(fill = "Covariate") +
  coord_polar("y")


for (i in c(2:6, 8:11)) {
  df_var <- read_csv(file.path(fp_evals, paste0("var_importance_", i, ".csv"))) %>%
    dplyr::select(contains('RUN'))
  
  rf_temp <- df_var %>% dplyr::select(contains("RF"))
  rf_temp <- rf_temp %>% gather() %>% mutate(Var = rep(env_covars,10),
                                                 Month = i)
  names(rf_temp) <- c("Model", "Perc.Cont", "Var", "Month") 
  rf_temp <- rf_temp %>% mutate(Model = "RF", 
                                    Perc.Cont.Temp = rep(rowMeans(df_var %>% 
                                                                    dplyr::select(
                                                                      contains("RF")), 
                                                                      na.rm = TRUE), 10)) %>%
    mutate(Perc.Cont.Norm = Perc.Cont/sum(unique(Perc.Cont.Temp))/10)
  
  new_df_rf <- rbind(new_df_rf, rf_temp)
  
  print(sum(rf_temp$Perc.Cont.Norm))
}

#Color blind friendly palette. Source: https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/#use-a-colorblind-friendly-palette
cbp <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(data = new_df_rf, aes(x = Month, fill = factor(Var, env_covars))) +
  geom_col(aes(y = Perc.Cont.Norm)) +
  scale_x_continuous(breaks = seq(1,12,2), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) + 
  scale_fill_viridis_d() +
  labs(x = "Month", y = "Contribution", fill = "Covariate") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))

avg <- new_df_rf %>%
  group_by(Var) %>%
  summarize(mean = mean(Perc.Cont)) %>%
  mutate(norm = mean/sum(mean))

