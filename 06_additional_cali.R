setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# clean up
rm(list=ls())
graphics.off()
cat("\14")

library(tidyverse)
library(grid)
library(gridExtra)
library(ggExtra)
library(ggdendro)
library(cluster)
library(ggplot2)
library(ggpubr)
library(knitr)
library(data.table)
library(viridis)

# source settings script
source("0_Settings.R")
# the results from the calibration runs descriptions of the different columns
# can be found in "data/results_lhc_description.csv"
res_all <- read.csv("data/results_lhc.csv")
# meta information about the lakes 
lake_meta <- readRDS("data_derived/lake_meta_data_derived.RDS")
lake_meta_desc <- readRDS("data_derived/lake_meta_desc_derived.RDS")

# data frame with all metrics for single best model per lake
best_all_a <- readRDS("data_derived/single_best_model.RDS")

# data frame with all metrics for best set per lake and per model
best_all <- readRDS("data_derived/best_par_sets.RDS")

# load additional calibration results
res_cali <- read.csv("data/results_lhc_revision.csv")

##------------- pick best performing parameter sets for the re run models -------------------------

# select best performing parameter set for each performance metric, lake,
# and model
# extract best (rmse) parameter set for each lake and model
best_rmse2 <- res_cali |> group_by(lake = lake,
                             model = model) |>
  slice_min(rmse) |> mutate(best_met = "rmse")
# extract best (r) parameter set for each lake and model
best_r2 <- res_cali |> group_by(lake = lake,
                          model = model) |>
  slice_max(r) |> mutate(best_met = "r")

# extract best (nse) parameter set for each lake and model
best_nse2 <- res_cali |> group_by(lake = lake,
                            model = model) |>
  slice_max(nse) |> mutate(best_met = "nse")

# extract best (bias) parameter set for each lake and model
best_bias2 <- res_cali |> group_by(lake = lake,
                             model = model) |>
  slice_min(abs(bias)) |> mutate(best_met = "bias")

# extract best (nmae) parameter set for each lake and model
best_nmae2 <- res_cali |> group_by(lake = lake,
                             model = model) |>
  mutate(nmae = ifelse(nmae > 1e5, 0, nmae)) |>
  filter(nmae != 0 & !is.infinite(nmae)) |>
  slice_min(nmae) |> mutate(best_met = "nmae")

# extract best (mae) parameter set for each lake and model
best_mae2 <- res_cali |> group_by(lake = lake,
                            model = model) |>
  slice_min(mae) |> mutate(best_met = "mae")

# data frame with all metrics for best set per lake, model, and metric
best_all2 <- rbind(best_bias2, best_mae2, best_nmae2,
                  best_nse2, best_r2, best_rmse2) |> ungroup()
rm(best_bias2, best_mae2, best_nmae2,
   best_nse2, best_r2, best_rmse2)

# add columns with missing cols
mcols <- colnames(best_all)[!(colnames(best_all) %in% colnames(best_all2))]
best_all2[mcols] <- NA

# filter only the performance metrics which we set in 0_Settings.R
best_all2 <- filter(best_all2, best_met %in% p_metrics) |>
  select(-which(colnames(best_all2) %in% setdiff(c("bias",
                                                  "mae",
                                                  "nmae",
                                                  "nse",
                                                  "r",
                                                  "rmse"), p_metrics)))

# combine the new data with the old calibration runs
best_all <- rbind(best_all, best_all2) |> ungroup()


## filter new best_all to find single best model per lake
# extract best (rmse) parameter set for each lake
s_best_rmse2 <- best_all |> group_by(lake = lake) |>
  slice_min(rmse) |> mutate(best_met = "rmse")
# extract best (r) parameter set for each lake
s_best_r2 <- best_all |> group_by(lake = lake) |>
  slice_max(r) |> mutate(best_met = "r")

# extract best (nse) parameter set for each lake
s_best_nse2 <- best_all |> group_by(lake = lake) |>
  slice_max(nse) |> mutate(best_met = "nse")

# extract best (bias) parameter set for each lake
s_best_bias2 <- best_all |> group_by(lake = lake) |>
  slice_min(abs(bias)) |> mutate(best_met = "bias")



# data frame with all metrics for single best model per lake
best_all_a_r2 <- rbind(s_best_bias2,
                     s_best_nse2, s_best_r2, s_best_rmse2) |> ungroup() 
rm(s_best_bias2,
   s_best_nse2, s_best_r2, s_best_rmse2)

# filter only the performance metrics which we set in 0_Settings.R
best_all_a_r2 <- filter(best_all_a_r2, best_met %in% p_metrics) |>
  select(-which(colnames(best_all_a_r2) %in% setdiff(c("bias",
                                                     "mae",
                                                     "nmae",
                                                     "nse",
                                                     "r",
                                                     "rmse"), p_metrics))) |>
  distinct()


##--------------- plots -----------------------------------------------------------

best_all |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  filter(model %in% c("GOTM", "GOTM_r",
                      "Simstrat", "Simstrat_r")) |>
  select(lake, model, kmcluster, !!p_metrics, best_met) |>
  pivot_longer(!!p_metrics) |>
  mutate(best_met = ifelse(best_met == "bias", best_met, toupper(best_met)),
         name = ifelse(name == "bias", name, toupper(name))) |>
  slice(which(best_met == name)) |>
  ggplot() + geom_violin(aes(y = value, x = model,
                             fill = model)) +
  geom_jitter(aes(y = value,  x = model), height = 0,
              width = 0.125, size = 2.5, col = "grey42", alpha = 0.5) +
  thm + scale_fill_viridis_d("Model", option = "E", end = 0.9) +
  facet_grid(best_met~kmcluster, scales = "free_y") + theme(legend.position = "top") +
  xlab("Model") + ylab("") +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

# wind speed scaling
best_all |> left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(best_met = ifelse(best_met == "bias", best_met, toupper(best_met))) |>
  ggplot() + geom_hline(aes(yintercept = 1), lwd = 1.25, lty = "dashed",
                        col = "grey42") +
  geom_boxplot(aes(x = as.numeric(kmcluster), y = wind_speed,
                   fill = kmcluster)) +
  facet_grid(best_met~model) + scale_fill_viridis_d("Cluster") + thm +
  xlab("Cluster") + ylab("Calibrated wind scaling (-)") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave("Plots/dist_wind_scaling_cluster_model_revison.pdf", width = 16,
       height = 9, device = cairo_pdf)



# swr scaling
best_all |> left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(best_met = ifelse(best_met == "bias", best_met, toupper(best_met))) |>
  ggplot() + geom_hline(aes(yintercept = 1), lwd = 1.25, lty = "dashed",
                        col = "grey42") +
  geom_boxplot(aes(x = as.numeric(kmcluster), y = swr, fill = kmcluster)) +
  facet_grid(best_met~model) + scale_fill_viridis_d("Cluster") + thm +
  xlab("Cluster") + ylab("Calibrated swr scaling (-)") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave("Plots/dist_swr_scaling_cluster_model_revision.pdf", width = 16,
       height = 9, device = cairo_pdf)

# Kw scaling
best_all |> left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(best_met = ifelse(best_met == "bias", best_met, toupper(best_met))) |>
  ggplot() + geom_hline(aes(yintercept = 1), lwd = 1.25, lty = "dashed",
                        col = "grey42") +
  geom_boxplot(aes(x = as.numeric(kmcluster), y = Kw, fill = kmcluster)) +
  facet_grid(best_met~model) + scale_fill_viridis_d("Cluster") + thm +
  xlab("Cluster") + ylab("Calibrated Kw scaling (-)") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# difference between the new models
dat_diff <- best_all |> filter(model %in% c("GOTM", "GOTM_r",
                                            "Simstrat", "Simstrat_r")) |>
  select(lake, model, best_met, !!p_metrics) |>
  pivot_longer(!!p_metrics, names_to = "metric") |>
  pivot_wider(names_from = model, values_from = value,
              id_cols = c(lake, best_met, metric)) |>
  mutate(Simstrat = Simstrat - Simstrat_r,
         GOTM = GOTM - GOTM_r) |> select(-GOTM_r, - Simstrat_r) |>
  pivot_longer(cols = c(Simstrat, GOTM), values_to = "difference")

dat_diff |> filter(metric == best_met) |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  ggplot() +
  geom_hline(aes(yintercept = 0), lwd = 1.3, col = "grey") +
  geom_violin(aes(x = name, y = difference, fill = name)) +
  geom_point(aes(x = name, y = difference), col = "grey42") +
  scale_fill_manual("Model", values = viridis_pal(option = "C",
                                                  end = 0.9)(4)[3:4]) + thm +
  facet_grid(metric~kmcluster, scales = "free") +
  ylab("Difference") + xlab("Model")

ggsave("Plots/diff_gotm_simstrat.pdf", width = 13, height = 7, device = cairo_pdf)

best_all |> filter(model %in% c("Simstrat", "GOTM")) |>
  pivot_longer(!!p_metrics, names_to = "metric") |>
  left_join(rename(dat_diff, "model" = "name"),
            by = c("lake", "model", "best_met", "metric")) |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  filter(metric == best_met) |>
  select(-par_id, -swr, -wind_speed, -Kw, - hgeo, - a_seiche, -c_relax_C) |>
  ggplot() + geom_point(aes(x = value, y = difference, col = model)) +
  facet_wrap(metric ~ kmcluster, scales = "free") +
  scale_color_manual("Model", values = viridis_pal(option = "C",
                                                  end = 0.9)(4)[3:4]) +
  thm


# distribution of performance for the models
best_all |> filter(model %in% c("GOTM", "GOTM_r",
                                "Simstrat", "Simstrat_r")) |>
  pivot_longer(!!p_metrics) |>
  filter(best_met == name) |>
  ggplot() +
  geom_histogram(aes(x = value, y = after_stat(count)),
                 bins = 30, col = 1) +
  thm + xlab("") +
  facet_grid(model~best_met, scales = "free")

# same plot but with panels for model combined and revised version in different
# color
best_all |> filter(model %in% c("GOTM", "GOTM_r",
                                "Simstrat", "Simstrat_r")) |>
  pivot_longer(!!p_metrics) |>
  filter(best_met == name) |>
  mutate(revised = grepl("\\_r", model),
         model = gsub("\\_r", "", model)) |>
  ggplot() +
  geom_histogram(aes(x = value, y = after_stat(count), fill = revised),
                 bins = 30, position = "dodge", col = 1) +
  thm + xlab("") +
  facet_grid(model~best_met, scales = "free")

##--------- plot with pie charts, but Simstrat with aseiche = 0

# calculate fraction of lakes for which each model performs best across the
# different metrics

count_best <- res_all |> select(lake, model, rmse, nse, r, bias, mae, nmae) |>
  rbind(select(res_cali, lake, model, rmse, nse, r, bias, mae, nmae)) |>
  filter(model %in% c("FLake", "GLM", "GOTM", "Simstrat_r")) |>
  group_by(lake) |>
  reframe(rmse = model[which.min(rmse)],
          nse = model[which.max(nse)],
          r = model[which.max(r)],
          bias = model[which.min(abs(bias))],
          mae = model[which.min(mae)],
          nmae = model[which.min(nmae)]) |>
  pivot_longer(-1) |> 
  filter(name %in% p_metrics) |>
  group_by(value, name) |>
  reframe(n = round(n()/73, 3)*100) |>
  rename(Metric = "name", Model = "value") |>
  pivot_wider(id_cols = "Model", names_from = "Metric", values_from = "n") |>
  setNames(c("Model", "bias", "mae", "nmae",
             "NSE", "r", "RMSE") [c(1, which(c("bias", "mae",
                                               "nmae", "nse",
                                               "r", "rmse") %in% p_metrics)+1)])



## filter new best_all to find single best model per lake
# extract best (rmse) parameter set for each lake
s_best_rmse2 <- best_all |>
  filter(model %in% c("FLake", "GLM", "GOTM", "Simstrat_r")) |>
  group_by(lake = lake) |>
  slice_min(rmse) |> mutate(best_met = "rmse")
# extract best (r) parameter set for each lake
s_best_r2 <- best_all |>
  filter(model %in% c("FLake", "GLM", "GOTM", "Simstrat_r")) |>
  group_by(lake = lake) |>
  slice_max(r) |> mutate(best_met = "r")

# extract best (nse) parameter set for each lake
s_best_nse2 <- best_all |>
  filter(model %in% c("FLake", "GLM", "GOTM", "Simstrat_r")) |>
  group_by(lake = lake) |>
  slice_max(nse) |> mutate(best_met = "nse")

# extract best (bias) parameter set for each lake
s_best_bias2 <- best_all |>
  filter(model %in% c("FLake", "GLM", "GOTM", "Simstrat_r")) |>
  group_by(lake = lake) |>
  slice_min(abs(bias)) |> mutate(best_met = "bias")



# data frame with all metrics for single best model per lake
best_all_as <- rbind(s_best_bias2,
                       s_best_nse2, s_best_r2, s_best_rmse2) |> ungroup() 
rm(s_best_bias2,
   s_best_nse2, s_best_r2, s_best_rmse2)

# filter only the performance metrics which we set in 0_Settings.R
best_all_as <- filter(best_all_as, best_met %in% p_metrics) |>
  select(-which(colnames(best_all_as) %in% setdiff(c("bias",
                                                       "mae",
                                                       "nmae",
                                                       "nse",
                                                       "r",
                                                       "rmse"), p_metrics))) |>
  distinct()

p <- list()
for(m in p_metrics) {
  
  p_dtmp <-  best_all_as |> pivot_longer(!!p_metrics) |>
    mutate(best_met = ifelse(best_met == "bias", best_met, toupper(best_met)),
           name = ifelse(name == "bias", name, toupper(name))) |>
    filter(best_met == name & best_met == ifelse(m == "bias", m, toupper(m))) |>
    ggplot() +
    geom_histogram(aes(x = value, y = after_stat(count)),
                   bins = 30, col = 1) +
    thm + xlab("") +
    facet_wrap(~best_met, scales = "free")
  
  p_pie <- count_best |> setNames(c("Model", p_metrics)) |>
    pivot_longer(!!m) |> arrange(rev(Model)) |> ggplot() +
    geom_col(aes(x = "", y = value, fill = Model),
             col = "white") +
    geom_text(aes(x = 1.8, y = value, label = paste0(value ,"%")),
              position = position_stack(vjust=0.46),
              size = 3.33, col = "black") +
    coord_polar("y", start = 0) + theme_void(base_size = 11) +
    labs(x = NULL, y = NULL, fill = NULL) +
    theme(legend.position = "right",
          plot.margin = margin(0,0,0,0),
          legend.key.height = unit(0.3, "cm"),
          legend.key.width = unit(0.3, "cm")) +
    scale_fill_viridis_d("Model", option = "C", end = 0.9)
  
  # combine the two previous plots
  
  xrng <- layer_scales(p_dtmp)$x$range$range
  yrng <- layer_scales(p_dtmp)$y$range$range
  
  p[[m]] <- p_dtmp + annotation_custom(ggplotGrob(p_pie),
                                       xmin = ifelse(m == "rmse",
                                                     mean(xrng) + 0.04*diff(xrng),
                                                     min(xrng) - 0.12*diff(xrng)),
                                       xmax = ifelse(m == "rmse",
                                                     max(xrng) + 0.12*diff(xrng),
                                                     mean(xrng) - 0.04*diff(xrng)),
                                       ymin = mean(yrng) - 0.175*diff(yrng),
                                       ymax = max(yrng) + 0.075*diff(yrng))
  
  
  
}

ggarrange(plotlist = p, ncol = 2, nrow = 2)


best_all |> left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  ggplot() + geom_point(aes(x = lake.area.sqkm, y = wind_speed, col = kmcluster)) +
  facet_grid(best_met~model) + scale_x_log10() + thm + scale_color_viridis_d("Cluster")



# relationship between k_min and lake size
best_all |> mutate(turb_param.k_min = ifelse(model == "GOTM_r", 1e-8, turb_param.k_min)) |>
  mutate(model = ifelse(model == "GOTM_r", "GOTM", model)) |>
  filter(model == "GOTM") |> group_by(lake) |> slice_min(rmse) |> filter(best_met == "rmse") |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  ggplot() + geom_point(aes(x = mean.depth.m*lake.area.sqkm, y = turb_param.k_min, col = kmcluster)) +
  thm + scale_y_log10() + scale_x_log10() +
  geom_smooth(method = "lm", aes(x = mean.depth.m*lake.area.sqkm, y = turb_param.k_min))

# relationship between a_seiche and lake size
best_all |> mutate(a_seiche = ifelse(model == "Simstrat_r", 0, a_seiche)) |>
  mutate(model = ifelse(model == "Simstrat_r", "Simstrat", model)) |>
  filter(model == "Simstrat") |> group_by(lake) |> slice_min(rmse) |> filter(best_met == "rmse") |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  ggplot() + geom_point(aes(x = mean.depth.m*lake.area.sqkm, y = a_seiche, col = kmcluster)) +
  thm + scale_y_log10() + scale_x_log10() +
  geom_smooth(method = "lm", aes(x = mean.depth.m*lake.area.sqkm, y = a_seiche))

