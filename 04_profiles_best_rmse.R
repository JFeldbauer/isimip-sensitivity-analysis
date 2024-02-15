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
  group_by(model, depth_bin, kmcluster) |>
  reframe(sd_rmse = sd(rmse, na.rm = TRUE),
          mean_rmse = mean(rmse, na.rm = TRUE),
          median_rmse = median(rmse, na.rm = TRUE),
          q25_rmse = quantile(rmse, 0.25, na.rm = TRUE),
          q75_rmse = quantile(rmse, 0.75, na.rm = TRUE)) |>
  ggplot() +
  # geom_ribbon(aes(x = depth_bin, ymin =mean_rmse - sd_rmse,
  #                ymax = mean_rmse + sd_rmse, fill = model)) +
  geom_line(aes(x = depth_bin, y = median_rmse, col = model)) +
  geom_point(aes(x = depth_bin, y = median_rmse, col = model), size = 4) +
  geom_errorbar(aes(x = depth_bin, ymin = q25_rmse,
                    ymax = q75_rmse, col = model),
                    width = .025) +
    # geom_smooth(aes(y = rel_resid, x = rel_depth,
  #                 col = model, fill = model)) +
  scale_x_reverse() + coord_flip() +
  facet_grid(.~kmcluster) + thm +
  scale_color_viridis_d("Model", option = "C") +
  #scale_fill_viridis_d("Model", option = "C", alpha = 0.6) +
  ylab("RMSE (°C)") + xlab("Relative depth (-)")

ggsave("Plots/profiles_best_rmse.png", p_prfl, width = 13, height = 9)
