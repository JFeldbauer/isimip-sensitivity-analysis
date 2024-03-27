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
  group_by(n) |> reframe(freq = sum(n)) |> ungroup() |>
  mutate(freq = freq/sum(freq)) |>
  ggplot() + geom_col(aes(x = n, y = freq)) + thm +
  xlab("Number of different 'best' model") + ylab("Proportion (-)")

# look ath how the other performance metrics are doing compared on
# which one you use to select parameters
p_cor_met <- list()
for(b in unique(best_all$best_met)) {
  p_cor_met[[b]] <- ggpubr::ggarrange(
  best_all |> filter(best_met == b) |>
    ggplot() + geom_point(aes_string(x = b, y = "rmse", col = "model")) +
    thm + scale_color_viridis_d("Model", option = "C"),
  best_all |> filter(best_met == b) |>
    ggplot() + geom_point(aes_string(x = b, y = "r", col = "model")) +
    thm + scale_color_viridis_d("Model", option = "C"),
  best_all |> filter(best_met == b) |>
    ggplot() + geom_point(aes_string(x = b, y = "nse", col = "model")) +
    thm + scale_color_viridis_d("Model", option = "C"),
  best_all |> filter(best_met == b) |>
    ggplot() + geom_point(aes_string(x = b, y = "bias", col = "model")) +
    thm + scale_color_viridis_d("Model", option = "C"),
  ncol = 4, common.legend = TRUE)
}
p_cor_met <- ggpubr::ggarrange(plotlist = p_cor_met, nrow = 4,
                               common.legend = TRUE)

### median and below-treshold values for model performance
df_best_a_rmse <- best_all_a |>
  filter(best_met == "rmse")
df_best_a_r <- best_all_a |>
  filter(best_met == "r")

median(df_best_a_rmse$rmse)
median(df_best_a_r$r)

# Calculate percentage below a certain threshold
rmse_threshold <- 2

df_best_a_rmse |> 
  filter(rmse < rmse_threshold) |> 
  nrow() /
  nrow(df_best_a_rmse) * 100

df_best <- best_all |> 
  filter(best_met == "rmse")
setDT(df_best)
rmse_table = df_best[, .(median_rmse = median(rmse),
            perc_under_thresh = sum(rmse < rmse_threshold) / .N * 100),
        by = model]
kable(rmse_table, format = "pipe")

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


# distribution of the single best model per lake for all metrics
p_dist_lake_a <- best_all_a |> pivot_longer(!!p_metrics) |>
  filter(best_met == name) |>
  ggplot() +
  geom_histogram(aes(x = value, y = after_stat(count)),
                 bins = 30, col = 1) +
  thm + xlab("") +
  facet_wrap(~best_met, scales = "free")

ggsave("Plots/best_fit.png", p_dist_lake_a, width = 13, height = 7)

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
    geom_text(aes(x = 1.8, y = value, label = paste0(value ,"%")),
              position = position_stack(vjust=0.46),
              size = 2.75, col = "black") +
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

p_pie <- ggarrange(plots = p, ncol = 2, nrow = 2)
ggsave("Plots/best_fit_pie.pdf", p_pie, width = 13, height = 7)

# distribution of the best parameter per model and lake for all metrics
p_dist_lake <- best_all |> pivot_longer(!!p_metrics) |>
  filter(best_met == name) |>
  ggplot() +
  geom_histogram(aes(x = value, y = after_stat(count)),
                 bins = 50, col = 1) +
  thm + xlab("") +
  facet_grid(best_met~model, scales = "free")

ggsave("Plots/best_fit_model.png", p_dist_lake, width = 13, height = 7)

# same plot but with clusters
# distribution of the single best model per lake for all 6 metrics
best_all_a |> pivot_longer(!!p_metrics) |> filter(best_met == name) |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |> 
  ggplot() +
  geom_histogram(aes(x = value, y = after_stat(count), fill = kmcluster),
                 bins = 20, col = 1) +
  thm + xlab("") +
  facet_wrap(~best_met, scales = "free") +
  scale_fill_viridis_d("Cluster") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

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
        legend.position = "top") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave("Plots/mapisimip.pdf", width = 11, height = 7, bg = "white")

## function that plots relating min RMSE to lake characteristics
plot_meta <- function(data, measure = "rmse") {
  p_meta1 <- ggplot(data) + geom_point(aes_string(x = "mean.depth.m", y = measure, col = "model")) +
    scale_x_log10() + ggtitle(paste("best", measure, "vs. mean lake depth")) + thm + xlab("Mean lake depth (m)") +
    ylab(measure)
  
  p_meta2 <- ggplot(data) + geom_point(aes_string(x = "lake.area.sqkm", y = measure, col = "model")) +
    scale_x_log10() + ggtitle(paste("best", measure, "vs. lake area")) + thm + xlab("Lake Area (km²)") +
    ylab(measure)
  
  p_meta3 <- ggplot(data) + geom_point(aes_string(x = "latitude.dec.deg", y = measure, col = "model")) +
    ggtitle(paste("best", measure, "vs.latitude")) + thm + xlab("Latitude (°N)") +
    ylab(measure)
  
  p_meta4 <- ggplot(data) + geom_point(aes_string(x = "longitude.dec.deg", y = measure, col = "model")) +
    ggtitle(paste("best", measure, "vs.longitude")) + thm + xlab("Longitude (°E)") +
    ylab(measure)
  
  
  p_meta6 <- ggplot(data) + geom_point(aes_string(x = "elevation.m", y = measure, col = "model")) +
    ggtitle(paste("best", measure, "vs. elevation")) + thm + xlab("Lake elevation (m asl)") +
    ylab(measure)
  
  p_meta7 <- ggplot(data) + geom_point(aes_string(x = "kw", y = measure, col = "model")) +
    ggtitle(paste("best", measure, "vs. Secchi disk depth")) + thm + xlab("Average Kw value (1/m)") +
    ylab(measure)
  
  
  p_meta8 <- ggplot(data) + geom_point(aes_string(x = "months_meas", y = measure, col = "model")) +
    ggtitle(paste("best", measure, "vs. no. months with obs")) + thm + xlab("Months with observations (-)") +
    ylab(measure)
  
  
  p_meta10 <- ggplot(data) + geom_point(aes_string(x = "reldepth_median", y = measure, col = "model")) +
    ggtitle(paste("best", measure, "vs. relative depth of most obs")) + thm + xlab("center of relative depth (-)") +
    ylab(measure)
  
  
  p1_meta <- ggpubr::ggarrange(p_meta1,
                       p_meta2,
                       p_meta3,
                       p_meta6,
                       p_meta7,
                       p_meta4,
                       ncol = 3, nrow = 2, common.legend = TRUE)
  
  p2_meta <- ggpubr::ggarrange(p_meta8,
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
res |> left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  pivot_longer(!!p_metrics) |> ggplot() +
  geom_violin(aes(x = model, y = value, fill = model)) + scale_y_log10() +
  facet_grid(name~kmcluster, scales = "free") + thm + scale_y_log10() +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

# check if the same models perform good or bad along the four models
p_mcomp1 <- best_all |> pivot_longer(!!p_metrics) |> 
  mutate(best_met = ifelse(best_met == "bias", best_met, toupper(best_met)),
         name = ifelse(name == "bias", name, toupper(name))) |>
  group_by(lake, best_met, name) |> filter(best_met == name) |>
  reframe(range = diff(range(value)),
          sd = sd(value,),
          min = min(value),
          max = max(value)) |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  ggplot() + geom_boxplot(aes(y = sd, x = as.numeric(kmcluster),
                             fill = kmcluster)) +
  geom_point(aes(y = sd, x = as.numeric(kmcluster)),
             col = alpha("grey42", 0.666),
             size = 3.5) +
  facet_wrap(~name, scales = "free_y") + thm +
  scale_fill_viridis_d("Cluster") + ylab("standard deviation model performacne") +
  xlab("Cluster") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

p_mcomp2 <- best_all |> pivot_longer(!!p_metrics) |>
  mutate(best_met = ifelse(best_met == "bias", best_met, toupper(best_met)),
         name = ifelse(name == "bias", name, toupper(name))) |>
  group_by(lake, best_met, name) |> filter(best_met == name) |>
  reframe(range = diff(range(value)),
          sd = sd(value),
          best = case_when(name == "bias" ~ min(abs(value)),
                            name %in% c("mae", "nmae", "rmse") ~ min(value),
                            name %in% c("r", "nse") ~ max(value)),
          worst = case_when(name == "bias" ~ max(abs(value)),
                            name %in% c("mae", "nmae", "rmse") ~ max(value),
                            name %in% c("r", "nse") ~ min(value))) |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  ggplot() + geom_point(aes(x = worst, y = best, col = kmcluster),
                        size = 3) +
  geom_abline(aes(intercept = 0, slope = 1), col = 2, lty = 17) +
  facet_wrap(~name, scales = "free") + thm +
  scale_color_viridis_d("Cluster") + xlab("Poorest model performance") +
  ylab("Best model performance") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))


p_mcomp3 <- best_all |> pivot_longer(!!p_metrics) |>
  mutate(best_met = ifelse(best_met == "bias", best_met, toupper(best_met)),
         name = ifelse(name == "bias", name, toupper(name))) |>
  group_by(lake, best_met, name) |> filter(best_met == name) |>
  reframe(range = diff(range(value)),
          sd = sd(value),
          best = case_when(name == "bias" ~ min(abs(value)),
                           name %in% c("mae", "nmae", "rmse") ~ min(value),
                           name %in% c("r", "nse") ~ max(value)),
          worst = case_when(name == "bias" ~ max(abs(value)),
                            name %in% c("mae", "nmae", "rmse") ~ max(value),
                            name %in% c("r", "nse") ~ min(value))) |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  ggplot() + geom_point(aes(x = best, y = sd, col = kmcluster),
                        size = 3) +
  facet_wrap(~name, scales = "free") + thm +
  scale_color_viridis_d("Cluster") + ylab("standard deviation model performacne") +
  xlab("Best model performance") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))


ggsave("Plots/range_best_model.pdf", p_mcomp1, width = 13, height = 9,
       bg = "white")

ggsave("Plots/poorest_best_model.png", p_mcomp2, width = 13, height = 9,
       bg = "white")

ggsave("Plots/best_model_range.png", p_mcomp3, width = 13, height = 9,
       bg = "white")

## just look at the deep and medium lakes to see which models perform better

p_dl <- best_all |> pivot_longer(!!p_metrics) |> 
  group_by(lake, best_met, name) |> filter(best_met == name) |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  filter(as.numeric(kmcluster) %in% 1:2) |>
  ggplot() + geom_boxplot(aes(y = value, x = model, fill = model)) +
    facet_wrap(~name, scales = "free") +
  scale_fill_viridis_d("Model", option = "C", end = 0.9) + thm +
  ggtitle("Deep and medium temperate lakes")

# look at the other cluster
p_ol <- best_all |> pivot_longer(!!p_metrics) |> 
  group_by(lake, best_met, name) |> filter(best_met == name) |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  filter(as.numeric(kmcluster) %in% 3:5) |>
  ggplot() + geom_boxplot(aes(y = value, x = model, fill = model)) +
  facet_wrap(~name, scales = "free") +
  scale_fill_viridis_d("Model", option = "C", end = 0.9) + thm +
  ggtitle("shallow, large shallow and warm lakes")

ggpubr::ggarrange(p_dl, p_ol, ncol = 2, common.legend = TRUE)
ggsave("Plots/deeper_shallower_perf.png", width = 13, height = 7)

##-------- relate calibrated parameter values to lake characteristics ----

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

ggsave("Plots/dist_wind_scaling_cluster_model.png", width = 13, height = 9)

best_all |> left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(best_met = ifelse(best_met == "bias", best_met, toupper(best_met))) |>
  ggplot() + geom_hline(aes(yintercept = 1), lwd = 1.25, lty = "dashed",
                        col = "grey42") +
  geom_boxplot(aes(x = model, y = wind_speed,
                   fill = model)) +
  facet_grid(best_met~kmcluster) +
  scale_fill_viridis_d("Model", option = "C", end = 0.9) + thm +
  xlab("Cluster") + ylab("Calibrated wind scaling (-)") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

ggsave("Plots/dist_wind_scaling_cluster.pdf", width = 13, height = 9)

# swr scaling
best_all |> left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(best_met = ifelse(best_met == "bias", best_met, toupper(best_met))) |>
  ggplot() + geom_hline(aes(yintercept = 1), lwd = 1.25, lty = "dashed",
                        col = "grey42") +
  geom_boxplot(aes(x = as.numeric(kmcluster), y = swr, fill = kmcluster)) +
  facet_grid(best_met~model) + scale_fill_viridis_d("Cluster") + thm +
  xlab("Cluster") + ylab("Calibrated swr scaling (-)") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave("Plots/dist_swr_scaling_cluster_model.png", width = 13, height = 9)

best_all |> left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(best_met = ifelse(best_met == "bias", best_met, toupper(best_met))) |>
  ggplot() + geom_hline(aes(yintercept = 1), lwd = 1.25, lty = "dashed",
                        col = "grey42") +
  geom_boxplot(aes(x = model, y = swr, fill = model)) +
  facet_grid(best_met~kmcluster) +
  scale_fill_viridis_d("Model", option = "C", end = 0.9) + thm +
  xlab("Cluster") + ylab("Calibrated swr scaling (-)") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

ggsave("Plots/dist_swr_scaling_cluster.png", width = 13, height = 9)

# kw scaling
best_all |> left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(best_met = ifelse(best_met == "bias", best_met, toupper(best_met))) |>
  ggplot() +
  geom_boxplot(aes(x = as.numeric(kmcluster), y = Kw, fill = kmcluster)) +
  facet_grid(best_met~model) +
  scale_fill_viridis_d("Cluster") + thm +
  xlab("Model") + ylab("Calibrated Kw scaling (-)") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_log10() 

ggsave("Plots/dist_kw_scaling_cluster_model.png", width = 13, height = 9)

best_all |> left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(best_met = ifelse(best_met == "bias", best_met, toupper(best_met))) |>
  ggplot() +
  geom_boxplot(aes(x = model, y = Kw, fill = model)) +
  facet_grid(best_met~kmcluster) +
  scale_fill_viridis_d("Model", option = "C", end = 0.9) + thm +
  xlab("Model") + ylab("Calibrated Kw scaling (-)") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

ggsave("Plots/dist_kw_scaling_cluster.png", width = 13, height = 9)

# all model specific parameter together
lapply(c("FLake", "GLM", "GOTM", "Simstrat"), function(m){
  best_all |> left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
    mutate(best_met = ifelse(best_met == "bias", best_met, toupper(best_met))) |>
    pivot_longer(cols = !!par_names) |> filter(model == m) |>
    filter(!(name %in% c("wind_speed", "swr", "Kw")))|>
    select(best_met, value, kmcluster, name) |> na.omit() |>
    ggplot() +
    geom_boxplot(aes(x = best_met, y = value, fill = kmcluster)) +
    facet_wrap(~name, scales = "free") + scale_fill_viridis_d("Cluster") +
    thm + xlab("") + ylab("") +
    guides(fill = guide_legend(nrow = 5, byrow = TRUE)) +
    theme(legend.position = "right")}) |>
  ggpubr::ggarrange(plotlist = _, labels = c("FLake", "GLM", "GOTM", "Simstrat"),
            common.legend = TRUE, ncol = 1, legend = "right")

ggsave("Plots/par_value_cluster.pdf", width = 17, height = 13, bg = "white")


# Distribution of scaling factors - RMSE only
Kw_range <- res |> group_by(lake, model) |>
  reframe(Kw_mid= (max(Kw) + min(Kw))/2)

best_all |> subset(best_met == "rmse") |>
  left_join(lake_meta, by = c("lake" = "Lake.Short.Name")) |>
  left_join(Kw_range) |>
  group_by(lake, model) |>
  reframe(Kw = Kw/Kw_mid,
          wind_speed = wind_speed,
          swr = swr,
          kmcluster = kmcluster) |>
  ungroup() |>
  pivot_longer(cols = c(wind_speed, swr, Kw)) |>
  ggplot() + geom_hline(aes(yintercept = 1), lwd = 1.25, lty = "dashed",
                        col = "grey42") +
  geom_boxplot(aes(x = model, y = value, fill = model)) +
  facet_grid(name~kmcluster, scales = "free") +
  scale_fill_viridis_d("Model", option = "C", end = 0.9) + thm +
  xlab("Cluster") + ylab("Scaling (-)") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

ggsave("Plots/dist_scaling_cluster_model_RMSE_only.pdf", width = 13, height = 9)

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
  
  p_i = ggpubr::ggarrange(plotlist = plts2, common.legend = T, align = "hv")
  
  plts[[length(plts) + 1]] = p_i
}
