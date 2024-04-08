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
library(lubridate)

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


# calculate ensemble mean
dat_ens <-  dat |> select(1:8) |>
  pivot_wider(id_cols = c(datetime, depth, lake, obs_temp),
              names_from = model, values_from = sim_temp) |>
  group_by(datetime, depth, lake) |>
  reframe(obs_temp = obs_temp,
          GLM = GLM,
          GOTM = GOTM,
          Simstrat = Simstrat,
          FLake = FLake,
          ensemble_mean = mean(c(GLM, GOTM, Simstrat, FLake), na.rm = TRUE)) |>
  pivot_longer(5:9, names_to = "model", values_to = "sim_temp") |>
  left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  left_join(max_dept) |>
  mutate(resid = sim_temp - obs_temp) |> mutate(rel_depth = depth/max_depth) |>
  mutate(rel_depth = ifelse(is.na(rel_depth), 0, rel_depth))


p_prfl <- dat |>
  mutate(depth_bin = cut(rel_depth,
                         breaks = seq(0, 1, 0.1),
                         labels = seq(0.05, 0.95, 0.1),
                         include.lowest = TRUE)) |>
  mutate(depth_bin = as.numeric(as.character(depth_bin))) |>
  group_by(model, lake, depth_bin, kmcluster, datetime) |>
  reframe(rmse = sqrt(mean(resid^2, na.rm = TRUE))) |>
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
                    ymin = q25_rmse,
                    ymax = q75_rmse, col = model),
                    width = .025, linewidth = 1.15) +
    # geom_smooth(aes(y = rel_resid, x = rel_depth,
  #                 col = model, fill = model)) +
  scale_x_reverse() + coord_flip() +
  facet_grid(.~kmcluster) + thm +
  scale_color_viridis_d("Model", option = "C", end = 0.9) +
  #scale_fill_viridis_d("Model", option = "C", alpha = 0.6) +
  ylab("RMSE (°C)") + xlab("Relative depth (-)")

ggsave("Plots/profiles_best_rmse.png", p_prfl, width = 13, height = 9)

  

## choose which function to use for thermocline depth
 thermo_fun <- function(...) thermo.depth(...)
# thermo_fun <- function(...) center.buoyancy(...)

# calculate thermocline depth and N2 for each lake and observation
thrm_dpth <- dat |> 
  distinct(lake, datetime, depth, obs_temp) |> arrange(lake, datetime, depth) |>
  group_by(lake, datetime) |> filter(length(unique(depth)) > 2) |>
  reframe(thermo_depth = thermo_fun(obs_temp, depth),
          N2 = max(buoyancy.freq(obs_temp, depth), na.rm = TRUE))

# calculate RMSE for water temperature at the thermocline for each lake and model
dat_thermo <- dat |> left_join(thrm_dpth) |> group_by(lake, datetime) |>
  reframe(dist_thermo = depth - thermo_depth,
          obs_temp = obs_temp,
          sim_temp = sim_temp,
          model = model) |> filter(!is.na(dist_thermo)) |>
  group_by(lake, datetime) |> slice_min(abs(dist_thermo)) |>
  group_by(lake, model) |>
  reframe(rmse = sqrt(mean((obs_temp - sim_temp)^2, na.rm = TRUE))) |>
  left_join(meta, by = c("lake" = "Lake.Short.Name")) 

p_thermo <- dat_thermo |> ggplot() +
  geom_boxplot(aes(x = model, y = rmse, fill = model)) +
  facet_grid(.~kmcluster) +
  scale_fill_viridis_d("Model", option = "C", end = 0.9) + thm +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  ylab("RMSE (K)") + xlab("Model")

ggsave("Plots/rmse_at_thermocline.png", p_thermo, width = 11, height = 6)


# same plot but split for lakes with large N2
dat |> left_join(thrm_dpth) |> group_by(lake, datetime, model) |>
  reframe(rmse = sqrt(mean((obs_temp - sim_temp)^2, na.rm = TRUE)),
          N2 = N2,
          lN = ifelse(N2 < 1e-4, "N^2 < 1e-4", "N^2 >= 1e-4")) |>
  ggplot() +
  geom_boxplot(aes(x = model, y = rmse, fill = model)) +
  facet_grid(.~lN) +
  scale_fill_viridis_d("Model", option = "C", end = 0.9) + thm +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  ylab("RMSE (K)") + xlab("Model")

# calculate RMSE for thermocline depth (in m)
rmse_thermo <- dat |> 
  group_by(lake, datetime, model) |>
  filter(length(unique(depth)) > 2) |>
  reframe(thermo_depth_sim = thermo_fun(sim_temp, depth)) |>
  left_join(thrm_dpth) |> na.omit() |> group_by(lake, model) |>
  reframe(rmse_thrm = sqrt(mean((thermo_depth_sim - thermo_depth)^2))) |>
  left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  ggplot() + geom_boxplot(aes(x = model, y = rmse_thrm, fill = model)) +
  facet_grid(.~kmcluster) +
  scale_fill_viridis_d("Model", option = "C", end = 0.9) + thm +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  ylab("RMSE (m)") + xlab("Model") + scale_y_log10()

ggsave("Plots/rmse_thermocline.pdf", rmse_thermo, width = 11, height = 6)

p_both <- ggpubr::ggarrange(p_prfl + ggtitle("(A) normalized profiles of RMSE"),
                    p_thermo + xlab("") + ggtitle("(B) RMSE at thermocline depth"),
                    ncol = 1, common.legend = TRUE)

ggsave("Plots/profile_and_thermo_rmse.pdf", p_both, width = 13, height = 12)



## check if ensemble mean is good predictor
rmse_ens <- dat_ens |>
  group_by(lake, model) |>
  reframe(rmse = sqrt(mean((sim_temp - obs_temp)^2, na.rm = TRUE))) |>
  left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(model = ifelse(model == "ensemble_mean", "Ensemble mean", model))

count_best <- rmse_ens |> group_by(lake, kmcluster) |>
  reframe(rmse = model[which.min(rmse)]) |>
  pivot_longer(3) |> 
  group_by(value, name, kmcluster) |>
  reframe(n =n()) |> group_by(kmcluster) |>
  reframe(value = value, n = n/sum(n)) |>
  rename(Cluster = "kmcluster", Model = "value", Fraction = "n") 

p_cnt_ens <- count_best |>  arrange(rev(Model)) |> ggplot() +
  geom_col(aes(x = "", y = Fraction, fill = Model),
           col = "white") +
  geom_text(aes(x = 1.8, y = Fraction, label = paste0(round(Fraction*100, 1) ,"%")),
            position = position_stack(vjust=0.46),
            size = 2.75, col = "black") +
  coord_polar("y", start = 0) + theme_void(base_size = 11) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(legend.position = "right",
        plot.margin = margin(0,0,0,0),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm")) +
  facet_grid(.~Cluster) + ylab("") +
  scale_fill_manual("Model", values = c("black", viridis::plasma(4, , end = 0.9)))

p_viol_ens <- rmse_ens |> ggplot() + geom_violin(aes(x = model, y = rmse, fill = model), position = "dodge") +
  geom_jitter(aes(y = rmse,  x = model), height = 0,
              width = 0.125, size = 2.5, col = "grey42", alpha = 0.5) +
  facet_grid(.~kmcluster) + thm +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  scale_fill_manual("Model", values = c("black", viridis::plasma(4, , end = 0.9)))

ggpubr::ggarrange(p_cnt_ens, p_viol_ens, ncol = 1, common.legend = TRUE,
                  legend = "none")
