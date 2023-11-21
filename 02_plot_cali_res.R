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

# source settings script
source("0_Settings.R")
# the results from the calibration runs descriptions of the different columns
# can be found in "data/results_lhc_description.csv"
res <- read.csv("data/results_lhc.csv")
# meta information about the lakes 
lake_meta <- readRDS("data_derived/lake_meta_data_derived.RDS")
lake_meta_desc <- readRDS("data_derived/lake_meta_desc_derived.RDS")

# data frame with all metrics for single best model per lake
best_all_a <- readRDS("data_derived/single_best_model.RDS")

# data frame with all metrics for best set per lake and per model
best_all <- readRDS("data_derived/best_par_sets.RDS")


##---------- look at best performing model per lake and metric -----------------

# calculate fraction of lakes for which each model performs best across the
# different metrics

count_best <- res |> group_by(lake) |>
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

# create a table that can be copy and pasted into the quarto document
kable(count_best, format = "pipe")

# look at the distribution of number of different best performing models over the
# different metrics

res |> group_by(lake) |>
  reframe(rmse = model[which.min(rmse)],
          nse = model[which.max(nse)],
          r = model[which.max(r)],
          bias = model[which.min(abs(bias))],
          mae = model[which.min(mae)],
          nmae = model[which.min(nmae)]) |>
  pivot_longer(-1) |> filter(name %in% p_metrics) |>
  group_by(lake) |> reframe(n = length(unique(value))) |>
  ggplot() + geom_histogram(aes(x = n))



##--------------- statistical models for best model/performacne ----------------
# fit multinomial log-linear model 
dat <- best_all_a |> filter(best_met == "rmse") |>
  left_join(lake_meta, c("lake" = "Lake.Short.Name")) |>
  mutate(model = factor(model, model, model),
         Reservoir.or.lake. = factor(Reservoir.or.lake.,
         Reservoir.or.lake.,
         Reservoir.or.lake.)) |> ungroup()

dat$model <- relevel(dat$model, ref = "Simstrat")
test <- multinom(model ~ kw + elevation.m + max.depth.m +
                   lake.area.sqkm + latitude.dec.deg + longitude.dec.deg + vd,
                 data = dat)
summary(test)
exp(coef(test))

pdat <- expand_grid(kw = seq(min(dat$kw), max(dat$kw), length.out = 5),
                    #kw = median(dat$kw),
                    elevation.m = mean(dat$elevation.m),
                    #max.depth.m = seq(min(dat$max.depth.m),
                    #                  max(dat$max.depth.m), length.out = 20),
                    max.depth.m = median(dat$max.depth.m),
                    #lake.area.sqkm = mean(dat$lake.area.sqkm),
                    lake.area.sqkm = seq(min(dat$lake.area.sqkm),
                                         max(dat$lake.area.sqkm),
                                         length.out = 20),
                    latitude.dec.deg = mean(dat$latitude.dec.deg),
                    longitude.dec.deg = mean(dat$longitude.dec.deg),
                    vd = mean(dat$vd))
  

res_m <- predict(test, newdata = pdat, "probs")
res_m <- cbind(res_m, pdat) |> pivot_longer(1:4)

res_m |> ggplot() + geom_line(aes(x = lake.area.sqkm, y = value, col = name)) +
  facet_grid(.~kw)

## estimate performance metric based upon lake characteristics

ggplot(dat) + geom_point(aes(x = kw, y = lake.area.sqkm, col = model)) +
  scale_y_log10() + scale_x_log10()


## model best rsme from lake propertiers

rmse_m <- best_all_a |> filter(best_met == "rmse") |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  lm(formula = rmse ~ (kw + elevation.m + max.depth.m +
               lake.area.sqkm + latitude.dec.deg + longitude.dec.deg +
               reldepth_median + months_meas + osgood + vd))

rmse_m_step <- step(rmse_m)

summary(rmse_m)
summary(rmse_m_step)

##--------------- plots looking at the best performing parameter set -------------


# distribution of the single best model per lake for all 6 metrics
p_dist_lake_a <- best_all_a |> pivot_longer(!!p_metrics) |>
  filter(best_met == name) |>
  ggplot() +
  geom_histogram(aes(x = value, y = ..count..),
                 bins = 20, col = 1) +
  thm + xlab("") +
  facet_wrap(~best_met, scales = "free")

ggsave("Plots/best_fit.png", p_dist_lake_a, width = 13, height = 7)

# distribution of the best parameter per model and lake for all 6 metrics
p_dist_lake <- best_all |> pivot_longer(!!p_metrics) |>
  filter(best_met == name) |>
  ggplot() +
  geom_histogram(aes(x = value, y = ..count..),
                 bins = 20, col = 1) +
  thm + xlab("") +
  facet_grid(model~best_met, scales = "free")

ggsave("Plots/best_fit_model.png", p_dist_lake, width = 13, height = 7)

# same plot but with clusters
# distribution of the single best model per lake for all 6 metrics
best_all_a |> pivot_longer(!!p_metrics) |> filter(best_met == name) |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |> 
  ggplot() +
  geom_histogram(aes(x = value, y = ..count.., fill = kmcluster),
                 bins = 20, col = 1) +
  thm + xlab("") +
  facet_wrap(~best_met, scales = "free") +
  scale_fill_viridis_d("Cluster")

ggsave("Plots/best_fit_clust.png", width = 13, height = 7)

 # a map of the lakes with the location color coded according to the best
# performing model
world <- map_data("world")
best_all |> filter(best_met == "rmse") |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  group_by(lake) |> slice(which.min(rmse)) |> ggplot() +
  geom_map(
    data = world, map = world, fill = alpha("#b36b00", 0.666),
    aes(map_id = region)
  ) +
  geom_point(aes(y = latitude.dec.deg, x = longitude.dec.deg, fill = kmcluster),
             size = 4, pch = 23, color = "white") + ylim(-90, 90) +
  xlim(-180, 180) + xlab("Longitude") + ylab("Latitude") +
  theme_minimal(base_size = 17) + scale_fill_viridis_d("Cluster") +
  theme(panel.background = element_rect(fill = '#9fbfdf'),
        panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
        panel.grid.minor = element_blank(),
        legend.position = "top")

ggsave("Plots/mapisimip.png", width = 11, height = 7, bg = "white")

## function that plots relating min RMSE to lake characteristics
plot_meta <- function(data, measure = "rmse") {
  p_meta1 <- ggplot(data) + geom_point(aes_string(x = "mean.depth.m", y = measure, col = "model")) +
    scale_x_log10() + ggtitle("min RMSE vs. mean lake depth") + thm + xlab("Mean lake depth (m)") +
    ylab("RMSE (K)")
  
  p_meta2 <- ggplot(data) + geom_point(aes_string(x = "lake.area.sqkm", y = measure, col = "model")) +
    scale_x_log10() + ggtitle("min RMSE vs. lake area") + thm + xlab("Lake Area (km²)") +
    ylab("RMSE (K)")
  
  p_meta3 <- ggplot(data) + geom_point(aes_string(x = "latitude.dec.deg", y = measure, col = "model")) +
    scale_x_log10() + ggtitle("min RMSE vs.latitude") + thm + xlab("Latitude (°N)") +
    ylab("RMSE (K)")
  
  p_meta4 <- ggplot(data) + geom_point(aes_string(x = "longitude.dec.deg", y = measure, col = "model")) +
    scale_x_log10() + ggtitle("min RMSE vs.longitude") + thm + xlab("Longitude (°E)") +
    ylab("RMSE (K)")
  
  
  p_meta6 <- ggplot(data) + geom_point(aes_string(x = "elevation.m", y = measure, col = "model")) +
    scale_x_log10() + ggtitle("min RMSE vs. elevation") + thm + xlab("Lake elevation (m asl)") +
    ylab("RMSE (K)")
  
  p_meta7 <- ggplot(data) + geom_point(aes_string(x = "kw", y = measure, col = "model")) +
    ggtitle("min RMSE vs. Secchi disk depth") + thm + xlab("Average Kw value (1/m)") +
    ylab("RMSE (K)")
  
  
  p_meta8 <- ggplot(data) + geom_point(aes_string(x = "months_meas", y = measure, col = "model")) +
    ggtitle("min RMSE vs. no. months with obs") + thm + xlab("Months with observations (-)") +
    ylab("RMSE (K)")
  
  
  p_meta10 <- ggplot(data) + geom_point(aes_string(x = "reldepth_median", y = measure, col = "model")) +
    ggtitle("min RMSE vs. relative depth of most obs") + thm + xlab("center of relative depth (-)") +
    ylab("RMSE (K)")
  
  
  p1_meta <- ggarrange(p_meta1,
                       p_meta2,
                       p_meta3,
                       p_meta6,
                       p_meta7,
                       p_meta4,
                       ncol = 3, nrow = 2, common.legend = TRUE)
  
  p2_meta <- ggarrange(p_meta8,
                       p_meta10,
                       ncol = 2, nrow = 2, common.legend = TRUE)
  return(list(p1_meta, p2_meta))
}


# plot meta data vs best r for each lake and model
best_all |> filter(best_met == "r") |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |> 
  plot_meta(measure = "r")

# same plot but just for the best model per lake
best_all_a |> filter(best_met == "rmse") |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  plot_meta()


##-------- compare models -----------------------

# distribution of metrics along all 2000 parameter sets and 73 lakes
res |> pivot_longer(!!p_metrics) |> ggplot() +
  geom_violin(aes(x = model, y = value, fill = model)) + scale_y_log10() +
  facet_wrap(~name, scales = "free_y") + thm + scale_y_log10()

# check if the same models perform good or bad along the four models
p_mcomp1 <- best_all |> pivot_longer(!!p_metrics) |> 
  group_by(lake, best_met, name) |> filter(best_met == name) |>
  reframe(range = diff(range(value)),
          sd = sd(value,),
          min = min(value),
          max = max(value)) |>
  mutate(range = ifelse(range > 1e3, NA, range)) |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  ggplot() + geom_boxplot(aes(y = range, x = kmcluster, fill = kmcluster)) +
  facet_wrap(~name, scales = "free_y") + thm +
  scale_fill_viridis_d("Cluster") + ylab("Range model performacne") +
  xlab("Cluster") + scale_y_log10()

p_mcomp2 <- best_all |> pivot_longer(!!p_metrics) |>
  group_by(lake, best_met, name) |> filter(best_met == name) |>
  reframe(range = diff(range(value)),
          sd = sd(value),
          best = case_when(name == "bias" ~ min(abs(value)),
                            name %in% c("mae", "nmae", "rmse") ~ min(value),
                            name %in% c("r", "nse") ~ max(value)),
          worst = case_when(name == "bias" ~ max(abs(value)),
                            name %in% c("mae", "nmae", "rmse") ~ max(value),
                            name %in% c("r", "nse") ~ min(value))) |>
  mutate(best = ifelse(best > 1e3, NA, best),
         worst = ifelse(worst > 1e3, NA, worst)) |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  ggplot() + geom_point(aes(x = worst, y = best, col = kmcluster),
                        size = 3) +
  geom_abline(aes(intercept = 0, slope = 1), col = 2, lty = 17) +
  facet_wrap(~name, scales = "free") + thm +
  scale_color_viridis_d("Cluster") + xlab("Poorest model performance") +
  ylab("Best model performance") + scale_x_log10() + scale_y_log10()


ggsave("Plots/range_best_model.png", p_mcomp1, width = 13, height = 9,
       bg = "white")

ggsave("Plots/poorest_best_model.png", p_mcomp2, width = 13, height = 9,
       bg = "white")

##-------- relate calibrated parameter values to lake characteristics ----

# plot relating used wind speed scaling factor to lake area
# for each model
best_all_a |> left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  ggplot() + geom_point(aes(y = wind_speed,
                            x = lake.area.sqkm,
                            col = model), size = 2) +
  thm + facet_wrap(.~best_met) + scale_x_log10()

best_all |> left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  ggplot() + geom_hline(aes(yintercept = 1), lwd = 1.25, lty = "dashed",
                        col = "grey42") +
  geom_boxplot(aes(x = kmcluster, y = wind_speed, fill = kmcluster)) +
  facet_grid(model~best_met) + scale_fill_viridis_d("Cluster") + thm +
  xlab("Cluster") + ylab("Calibrated wind scaling (-)")

ggsave("Plots/dist_wind_scaling_cluster.png", width = 13, height = 11)

best_all |> left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  ggplot() + geom_hline(aes(yintercept = 1), lwd = 1.25, lty = "dashed",
                        col = "grey42") +
  geom_boxplot(aes(x = kmcluster, y = swr, fill = kmcluster)) +
  facet_grid(model~best_met) + scale_fill_viridis_d("Cluster") + thm +
  xlab("Cluster") + ylab("Calibrated swr scaling (-)")

ggsave("Plots/dist_swr_scaling_cluster.png", width = 13, height = 11)


lapply(c("FLake", "GLM", "GOTM", "Simstrat"), function(m){
  best_all |> left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
    pivot_longer(cols = !!par_names) |> filter(model == m) |>
    select(best_met, value, kmcluster, name) |> na.omit() |>
    ggplot() +
    geom_boxplot(aes(x = best_met, y = value, fill = kmcluster)) +
    facet_wrap(~name, scales = "free") + scale_fill_viridis_d("Cluster") +
    thm + xlab("") + ylab("")}) |>
  ggarrange(plotlist = _, labels = c("FLake", "GLM", "GOTM", "Simstrat"),
            common.legend = TRUE)

ggsave("Plots/par_value_cluster.png", width = 30, height = 20, bg = "white")


##### Plot best parameter values vs lake characteristics
lm <- read.csv("data/Lake_meta.csv")
df_best_rmse = best_all |> filter(best_met == "rmse") |>
  left_join(lm, by = c("lake" = "Lake.Short.Name")) |>data.table()


cal_pars = par_names
lake_chars = names(df_best_rmse)[c(28, 31, 33, 34, 37)]

plts = list()
for(i in lake_chars){
  plts2 = list()
  for(j in cal_pars){
    p = ggplot(df_best_rmse) +
      geom_point(aes(.data[[i]], .data[[j]], colour = model))
    if(i %in% c("mean.depth.m", "lake.area.sqkm")){
      p = p +
        scale_x_continuous(trans = "log10")
    }
    p = p +
      thm
    
    plts2[[length(plts2) + 1]] = p
  }
  
  p_i = ggarrange(plotlist = plts2, common.legend = T, align = "hv")
  
  plts[[length(plts) + 1]] = p_i
}
