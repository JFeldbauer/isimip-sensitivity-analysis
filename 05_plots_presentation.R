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
library(nnet)
library(data.table)
library(egg)

##-------------- read in data ----------------------------------------------
# source settings script
source("0_Settings.R")
# load results of sensitivity analysis
res <- read.csv("sensitivity_analysis/res_sens.csv")

# filter metrics that are set in 0_Settings.R
res <- filter(res, var %in% p_metrics)

# save oiriginal results for plots
res_o <- res
# load results of latin hypercube calibration
res_cali <- read.csv("data/results_lhc.csv")
# load lake meta data
meta <- readRDS("data_derived/lake_meta_data_derived.RDS")
meta_desc <- readRDS("data_derived/lake_meta_desc_derived.RDS")

# data frame with all metrics for best set per lake and per model
best_all <- readRDS("data_derived/best_par_sets.RDS")

# set sensitivity values that are smaller than the dummy variable
# to zero
res <- group_by(res, lake, model, var) |>
  reframe(delta = ifelse((delta - delta_conf) <=
                           (delta[names == "dummy"] + delta_conf[names == "dummy"]),
                         0, delta),
          S1 = ifelse((S1 - S1_conf) <=
                        (S1[names == "dummy"] + S1_conf[names == "dummy"]),
                      0, S1),
          names = names,
          delta_conf = delta_conf,
          S1_conf = S1_conf)


# data frame with all metrics for single best model per lake
best_all_a <- readRDS("data_derived/single_best_model.RDS")

# data frame with all metrics for best set per lake and per model
best_all <- readRDS("data_derived/best_par_sets.RDS")

# load profiles of best rmse runs
dat <- readRDS("data/all_profile_rmse_cali.RDS")

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

# larger base seize
thm <- theme_pubr(base_size = 23) + grids()

##---------- plot best model --------------------------------------

count_best <- res_cali |> group_by(lake) |>
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

# plot with pie charts
p <- list()
for(m in p_metrics) {
  
  p_dtmp <-  best_all_a |> pivot_longer(!!p_metrics) |>
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
    geom_text(aes(x = 1.7, y = value, label = paste0(value ,"%")),
              position = position_stack(vjust=0.46),
              size = 3.33, col = "black") +
    coord_polar("y", start = 0) + theme_void(base_size = 11) +
    labs(x = NULL, y = NULL, fill = NULL) +
    theme(legend.position = "none",
          plot.margin = margin(0,0,0,0),
          legend.key.height = unit(0.3, "cm"),
          legend.key.width = unit(0.3, "cm")) +
    scale_fill_viridis_d("Model", option = "C", end = 0.9)
  
  # combine the two previous plots
  
  xrng <- layer_scales(p_dtmp)$x$range$range
  yrng <- layer_scales(p_dtmp)$y$range$range
  
  p[[m]] <- p_dtmp + annotation_custom(ggplotGrob(p_pie),
                                       xmin = ifelse(m == "rmse",
                                                     mean(xrng) + 0.01*diff(xrng),
                                                     min(xrng) - 0.22*diff(xrng)),
                                       xmax = ifelse(m == "rmse",
                                                     max(xrng) + 0.22*diff(xrng),
                                                     mean(xrng) - 0.01*diff(xrng)),
                                       ymin = mean(yrng) - 0.24*diff(yrng),
                                       ymax = max(yrng) + 0.14*diff(yrng))
  
  
  
}

p_leg <- count_best |> setNames(c("Model", p_metrics)) |>
  pivot_longer(!!m) |> arrange(rev(Model)) |> ggplot() +
  geom_col(aes(x = "", y = value, fill = Model),
           col = "white") + thm +
  theme(legend.position = "bottom",
        plot.margin = margin(0,0,0,0),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.6, "cm")) +
  scale_fill_viridis_d("Model", option = "C", end = 0.9)
p_leg <- get_legend(p_leg)
p_leg <- as_ggplot(p_leg)

p_pie <- ggarrange(plots = p, ncol = 2, nrow = 2)
p_pie <- ggarrange(as_ggplot(p_pie), p_leg, nrow = 2, heights = c(0.9, 0.1))
ggsave("Plots/Plots_DGL/best_fit_pie.pdf", p_pie, width = 13, height = 9)

##-------- plot map --------------------------------------------

world <- map_data("world")
best_all |> filter(best_met == "rmse") |>
  left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  group_by(lake) |> slice(which.min(rmse)) |> ggplot() +
  geom_map(
    data = world, map = world, fill = alpha("#b36b00", 0.666),
    aes(map_id = region)
  ) +
  geom_point(aes(y = latitude.dec.deg, x = longitude.dec.deg, fill = kmcluster),
             size = 4, pch = 23, color = "white") + ylim(-90, 90) +
  xlim(-180, 180) + xlab("Longitude") + ylab("Latitude") +
  theme_minimal(base_size = 26) + scale_fill_viridis_d("Cluster") +
  theme(panel.background = element_rect(fill = '#9fbfdf'),
        panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
        panel.grid.minor = element_blank(),
        legend.position = "top") +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

ggsave("Plots/Plots_DGL/mapisimip.pdf", width = 12, height = 9, bg = "white")

##---------- plot profiles -----------------------------------------------------

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
  facet_wrap(.~kmcluster) + thm +
  scale_color_viridis_d("Model", option = "C", end = 0.9) +
  #scale_fill_viridis_d("Model", option = "C", alpha = 0.6) +
  ylab("RMSE (K)") + xlab("Relative depth (-)")

ggsave("Plots/Plots_DGL/rmse_profile.pdf", plot = p_prfl, width = 17,
       height = 13, bg = "white")

##----------- ensemble mean ---------------------------------------------

## check if ensemble mean is good predictor
rmse_ens <- dat_ens |>
  group_by(lake, model) |>
  reframe(rmse = sqrt(mean((sim_temp - obs_temp)^2, na.rm = TRUE))) |>
  left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(model = ifelse(model == "ensemble_mean", "Ens. mean", model))

count_best <- rmse_ens |> group_by(lake, kmcluster) |>
  reframe(rmse = model[which.min(rmse)]) |>
  pivot_longer(3) |> 
  group_by(value, name, kmcluster) |>
  reframe(n =n()) |> group_by(kmcluster) |>
  reframe(value = value, frac = n/sum(n), n = n) |>
  rename(Cluster = "kmcluster", Model = "value", Fraction = "frac") 

p_cnt_ens <- count_best |>
  mutate(Model = factor(Model, levels = c("Ens. mean", "FLake", "GLM", "GOTM","Simstrat")),
         Mod_f = abs(as.numeric(Model) - 6)) |>
  arrange(Cluster, Mod_f) |>
  ggplot() +
  geom_col(aes(x = "", y = Fraction, fill = Model),
           col = "white") +
  geom_text(aes(x = 1.73, y = Fraction, label = paste0(round(Fraction*100, 1) ,"%")),
            position = position_stack(vjust=0.5),
            size = 5.5, col = "grey42") +
  coord_polar("y", start = 0) + theme_void(base_size = 23) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  labs(x = " ", y = "", fill = NULL) +
  theme(legend.position = "right",
        plot.margin = margin(0, 0.5, 0, 3.5, unit = "lines"),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip
        panel.grid = element_blank()) +
  facet_grid(.~Cluster) + ylab("") +
  scale_fill_manual("Model", values = c("black", viridis::plasma(4, end = 0.9)))

p_viol_ens <- rmse_ens |> ggplot() + geom_violin(aes(x = model, y = rmse, fill = model), position = "dodge") +
  geom_jitter(aes(y = rmse,  x = model), height = 0,
              width = 0.125, size = 2.5, col = "grey42", alpha = 0.5) +
  facet_grid(.~kmcluster) + thm +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1),
        plot.margin = margin(0, 0.5, 0, 0.5, unit = "lines")) +
  scale_fill_manual("Model", values = c("black", viridis::plasma(4, end = 0.9))) +
  xlab("Model") + ylab("RMSE (K)")

ggpubr::ggarrange(p_cnt_ens, p_viol_ens, ncol = 1, common.legend = TRUE,
                  legend = "none")
ggsave("Plots/Plots_DGL/count_best_ensemble_mean.pdf", width = 16, height = 10.5)

##----------------- plot performacnce per cluster ------------------------------

p_rmsec <- best_all_a |>
  left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  select(lake, model, kmcluster, !!p_metrics, best_met) |>
  pivot_longer(!!p_metrics) |>
  mutate(best_met = ifelse(best_met == "bias", best_met, toupper(best_met)),
         name = ifelse(name == "bias", name, toupper(name))) |>
  slice(which(best_met == name)) |>
  ggplot() + geom_violin(aes(y = value, x = as.numeric(kmcluster),
                             fill = kmcluster)) +
  geom_jitter(aes(y = value,  x = as.numeric(kmcluster)), height = 0,
              width = 0.125, size = 2.5, col = "grey42", alpha = 0.5) +
  thm + scale_fill_viridis_d("Cluster") +
  facet_wrap(~best_met, scales = "free_y") + theme(legend.position = "top") +
  xlab("Cluster") + ylab("") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave("Plots/Plots_DGL/performance_cluster.pdf", p_rmsec, width = 13, height = 9, bg = "white")


##-------------- plot cluster characteristics --------------------------------------

dat_clust <- meta |>
  select(-Lake.Name, -Lake.Name.Folder, -Country, -Average.Secchi.disk.depth.m,
         -Light.extinction.coefficient.m, 
         -depth_meas, -kw_sd, -tsurf_sd, -tbot_sd) |>
  #-tbot, -min_tsurf, -max.depth.m, -osgood, -months_median, -Reservoir.or.lake.) |>
  rename(lake = "Lake.Short.Name")



# distributuin of the lake characteristics
p_clst_char <- lapply(c(colnames(dat_clust)[!colnames(dat_clust) %in% c("lake",
                                                                        "kmcluster",
                                                                        "hcluster")],
                        "legend"),
                      function(c) {
                        if(c != "legend"){
                          dat <- select(dat_clust, c, "kmcluster") |> mutate(kmclustern = as.numeric(kmcluster))
                          desc <- meta_desc$short_description[meta_desc$column_name == c]
                          unit <- meta_desc$unit[meta_desc$column_name == c]
                          if(is.factor(dat[, c])) {
                            p <- dat |> table() |> as.data.frame() |> group_by(kmcluster) |>
                              ggplot() +
                              geom_col(aes_string(fill = c, y = "Freq", x = "kmclustern")) +
                              # scale_fill_viridis_d(desc, option = ifelse(c == "vd", "G", "E")) +
                              xlab("") + thm +
                              theme(legend.position = "top",
                                    legend.title=element_blank()) + guides(fill=guide_legend(ncol=2))
                          } else {
                            p <- dat |> ggplot() +
                              geom_violin(aes_string(y = c, fill = "kmcluster", x = "kmclustern")) +
                              geom_jitter(aes_string(y = c,  x = "kmclustern"), height = 0,
                                          width = 0.125, size = 2.5, col = "grey42", alpha = 0.5) +
                              scale_fill_viridis_d("") +
                              thm + xlab("") +
                              ylab(paste0(desc, " ( ", unit, " )")) +
                              theme(legend.position = "none")
                          }
                          
                          if(c %in% c("max.depth.m",
                                      "mean.depth.m",
                                      "lake.area.sqkm")) {
                            p <- p + scale_y_log10()
                          }
                        } else {
                          dat <- data.frame(x = 0,
                                            y = rev(1:5),
                                            leg = paste(1:5, "-",levels(meta$kmcluster)))
                          p <- dat |> ggplot() +
                            geom_point(aes(x = x, y = y, col = leg), size = 6.66) +
                            geom_text(aes(x = x + 0.1, y = y, label = leg),
                                      size = 4, hjust = 0) +
                            theme_void() + xlim(-0.15, 0.85) + ylim(0, 5) +
                            theme(legend.position = "none") + scale_color_viridis_d("")
                          
                        }
                        
                        return(p)
                      }) |> ggpubr::ggarrange(plotlist = _, ncol = 6, nrow = 3)

ggsave("Plots/Plots_DGL/clust_char.pdf", p_clst_char, width = 21, height = 13, 
       bg = "white", device = cairo_pdf)


##----------------- plots sensitivity --------------------------------------------

# boxplots of sensitivity metrics
res_o |> pivot_longer(cols = c(delta, S1)) |> filter(names != "dummy") |>
  mutate(name = case_match(name,
                           "delta" ~ "\u03B4",
                           "S1" ~ "S1")) |>
  mutate(names = gsub("mixing.coef_|turb_param.", "", names)) |>
  mutate(var = ifelse(var == "bias", var, toupper(var))) |>
  filter(var %in% c("RMSE", "R")) |>
  ggplot() +
  geom_boxplot(aes(x = names, y = value, fill = name)) +
  facet_grid(var ~ model, scales = "free") +
  scale_fill_manual("Sensitivity measure",
                    values = c("#45B2DD", "#72035F")) +
  thm + xlab("") +
  theme(axis.text.x=element_text(angle = -60, hjust = 0)) 

ggsave("Plots/Plots_DGL/sens_value.pdf", width = 12, height = 10, device = cairo_pdf)
