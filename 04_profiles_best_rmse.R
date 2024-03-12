setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# clean up
rm(list = ls())
graphics.off()
cat("\14")


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggExtra)
library(ggdendro)
library(rLakeAnalyzer)

##-------------- read in data ----------------------------------------------
# source settings script
source("0_Settings.R")
# load profiles of best rmse runs
dat <- readRDS("data/all_profile_rmse_cali.RDS")

# load lake meta data
meta <- readRDS("data_derived/lake_meta_data_derived.RDS")
meta_desc <- readRDS("data_derived/lake_meta_desc_derived.RDS")

# get maximum depth where there was a observation for each lake
max_dept <- dat |> group_by(lake) |> reframe(max_depth = max(depth))

dat <- left_join(dat, meta, by = c("lake" = "Lake.Short.Name")) |>
  left_join(max_dept) |>
  mutate(resid = sim_temp - obs_temp) |> mutate(rel_depth = depth/max_depth) |>
  mutate(rel_depth = ifelse(is.na(rel_depth), 0, rel_depth))


# test if the result is consistent
dat |> group_by(lake, model, kmcluster) |> reframe(rmse = sqrt(mean(resid^2))) |>
  ggplot() + geom_boxplot(aes(x = model, y = rmse, fill = model)) +
  facet_wrap(.~kmcluster) + thm +
  scale_fill_viridis_d("Model", option = "C")

p_prfl <- dat |> mutate(depth_bin = cut(rel_depth,
                                        breaks = seq(0, 1, 0.1),
                                        labels = seq(0.05, 0.95, 0.1),
                                        include.lowest = TRUE)) |>
  mutate(depth_bin = as.numeric(as.character(depth_bin))) |>
  group_by(model, lake, depth_bin, kmcluster, datetime) |>
  reframe(rmse = sqrt(mean(resid^2))) |>
  group_by(model, depth_bin, kmcluster, lake) |>
  reframe(rmse = median(rmse)) |>
  group_by(model, depth_bin, kmcluster) |>
  reframe(sd_rmse = sd(rmse, na.rm = TRUE),
          mean_rmse = mean(rmse, na.rm = TRUE),
          median_rmse = median(rmse, na.rm = TRUE),
          q25_rmse = quantile(rmse, 0.25, na.rm = TRUE),
          q75_rmse = quantile(rmse, 0.75, na.rm = TRUE),
          q05_rmse = quantile(rmse, 0.05, na.rm = TRUE),
          q95_rmse = quantile(rmse, 0.95, na.rm = TRUE),
          min_rmse = min(rmse, na.rm = TRUE),
          max_rmse = max(rmse, na.rm = TRUE)) |>
  ggplot() +
  # geom_ribbon(aes(x = depth_bin, ymin =mean_rmse - sd_rmse,
  #                ymax = mean_rmse + sd_rmse, fill = model)) +
  geom_line(aes(x = depth_bin + (as.numeric(as.factor(model))-1)/80,
                y = median_rmse, col = model), lty = "dashed") +
  geom_point(aes(x = depth_bin + (as.numeric(as.factor(model))-1)/80,
                 y = median_rmse, col = model), size = 4) +
  geom_errorbar(aes(x = depth_bin + (as.numeric(as.factor(model))-1)/80,
                    ymin = q05_rmse,
                    ymax = q95_rmse, col = model),
                    width = .025, linewidth = 1.15) +
    # geom_smooth(aes(y = rel_resid, x = rel_depth,
  #                 col = model, fill = model)) +
  scale_x_reverse() + coord_flip() +
  facet_grid(.~kmcluster) + thm +
  scale_color_viridis_d("Model", option = "C") +
  #scale_fill_viridis_d("Model", option = "C", alpha = 0.6) +
  ylab("RMSE (°C)") + xlab("Relative depth (-)")

ggsave("Plots/profiles_best_rmse.png", p_prfl, width = 13, height = 9)


# calculate thermocline depth for each lake and observation
thrm_dpth <- dat |> 
  distinct(lake, datetime, depth, obs_temp) |>
  group_by(lake, datetime) |>
  reframe(thermo_depth = thermo.depth(obs_temp, depth))

# calculate RMSE for water temperature at the thermocline for each lake and model
dat_thermo <- dat |> left_join(thrm_dpth) |> group_by(lake, datetime) |>
  reframe(dist_thermo = depth - thermo_depth,
          obs_temp = obs_temp,
          sim_temp = sim_temp,
          model = model) |> filter(!is.na(dist_thermo)) |>
  group_by(lake, datetime) |> slice_min(abs(dist_thermo)) |>
  group_by(lake, model) |>
  reframe(rmse = sqrt(mean((obs_temp - sim_temp)^2))) |>
  left_join(meta, by = c("lake" = "Lake.Short.Name"))

p_thermo <- dat_thermo |> ggplot() + geom_boxplot(aes(x = model, y = rmse, fill = model)) +
  facet_grid(.~kmcluster) + scale_fill_viridis_d("Model", "C") + thm +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  ylab("RMSE (K)") + xlab("Model")

ggsave("Plots/rmse_at_thermocline.png", p_thermo, width = 11, height = 6)

# calculate RMSE for thermocline depth (in m)
rmse_thermo <- dat |> 
  group_by(lake, datetime, model) |>
  reframe(thermo_depth_sim = thermo.depth(sim_temp, depth)) |>
  left_join(thrm_dpth) |> na.omit() |> group_by(lake, model) |>
  reframe(rmse_thrm = sqrt(mean((thermo_depth_sim - thermo_depth)^2))) |>
  left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  ggplot() + geom_boxplot(aes(x = model, y = rmse_thrm, fill = model)) +
  facet_grid(.~kmcluster) + scale_fill_viridis_d("Model", "C") + thm +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  ylab("RMSE (m)") + xlab("Model")

ggsave("Plots/rmse_thermocline.png", rmse_thermo, width = 11, height = 6)
