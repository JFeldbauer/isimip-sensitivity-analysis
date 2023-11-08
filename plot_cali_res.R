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

# the results from the calibration runs descriptions of the different columns
# can be found in "data/results_lhc_description.csv"
res <- read.csv("data/results_lhc.csv")
# meta information about the lakes a description of the columns can be
# found in "data/Lake_meta_description.csv"
lake_meta <- read.csv("data/Lake_meta.csv")


# ggplot theme to be used
thm <- theme_pubr(base_size = 15) + grids()


# extract best (rmse) parameter set for each lake and model
best_rmse <- res |> group_by(lake = lake,
                             model = model) |>
  slice_min(rmse) |> mutate(best_met = "rmse")
# extract best (r) parameter set for each lake and model
best_r <- res |> group_by(lake = lake,
                          model = model) |>
  slice_max(r) |> mutate(best_met = "r")

# extract best (nse) parameter set for each lake and model
best_nse <- res |> group_by(lake = lake,
                            model = model) |>
  slice_max(nse) |> mutate(best_met = "nse")

# extract best (bias) parameter set for each lake and model
best_bias <- res |> group_by(lake = lake,
                             model = model) |>
  slice_min(abs(bias)) |> mutate(best_met = "bias")

# extract best (nmae) parameter set for each lake and model
best_nmae <- res |> group_by(lake = lake,
                             model = model) |>
  slice_min(nmae) |> mutate(best_met = "nmae")

# extract best (mae) parameter set for each lake and model
best_mae <- res |> group_by(lake = lake,
                            model = model) |>
  slice_min(mae) |> mutate(best_met = "mae")


# extract best (rmse) parameter set for each lake
best_rmse_a <- res |> group_by(lake = lake) |>
  slice_min(rmse) |> mutate(best_met = "rmse")
# extract best (r) parameter set for each lake
best_r_a <- res |> group_by(lake = lake) |>
  slice_max(r) |> mutate(best_met = "r")

# extract best (nse) parameter set for each lake
best_nse_a <- res |> group_by(lake = lake) |>
  slice_max(nse) |> mutate(best_met = "nse")

# extract best (bias) parameter set for each lake
best_bias_a <- res |> group_by(lake = lake) |>
  slice_min(abs(bias)) |> mutate(best_met = "bias")

# extract best (nmae) parameter set for each lake
best_nmae_a <- res |> group_by(lake = lake) |>
  slice_min(nmae) |> mutate(best_met = "nmae")

# extract best (mae) parameter set for each lake
best_mae_a <- res |> group_by(lake = lake) |>
  slice_min(mae) |> mutate(best_met = "mae")

# data frame with all metrics for single best model per lake
best_all_a <- rbind(best_bias_a, best_mae_a, best_nmae_a,
                    best_nse_a, best_r_a, best_rmse_a)

# data frame with all metrics for best set per lake and per model
best_all <- rbind(best_bias, best_mae, best_nmae,
                  best_nse, best_r, best_rmse)


##---------- look at best performing model per lake and metric -----------------

# calculate fraction of lakes for which each model performs best across the
# different metrics

count_best <- res |> group_by(lake) |>
  reframe(rmse = model[which.min(rmse)],
          nse = model[which.max(nse)],
          r = model[which.max(r)],
          bias = model[which.max(abs(bias))],
          mae = model[which.min(mae)],
          nmae = model[which.min(nmae)]) |>
  pivot_longer(-1) |> group_by(value, name) |> reframe(n = round(n()/73, 3)*100) |>
  rename(Metric = "name", Model = "value") |>
  pivot_wider(id_cols = "Model", names_from = "Metric", values_from = "n") |>
  setNames(c("Model", "bias", "mae", "nmae",
             "NSE", "r", "RMSE"))

# create a table that can be copy and pasted into the quarto document
kable(count_best, format = "pipe")

# look at the distribution of number of different best performing models over the
# different metrics

res |> group_by(lake) |>
  reframe(rmse = model[which.min(rmse)],
          nse = model[which.max(nse)],
          r = model[which.max(r)],
          bias = model[which.max(abs(bias))],
          mae = model[which.min(mae)],
          nmae = model[which.min(nmae)]) |>
  pivot_longer(-1) |> group_by(lake) |> reframe(n = length(unique(value))) |>
  ggplot() + geom_histogram(aes(x = n))


##---------------- cluster analysis --------------------------------------------

# cluster analysis to try to estimate best performing model based on lake
# characteristics

# get average Kw values from all best measures for each lake
kw <- best_all |> group_by(lake) |> reframe(kw = mean(Kw))

# z-score normalize data for single best model
best_norm_a <- left_join(best_all_a, lake_meta,
                         by = c("lake" = "Lake.Short.Name")) |>
  left_join(kw) |>
  ungroup() |>
  mutate(across(c(4:24, 30:43),
                function(x)(x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))

# z-score normalize data for best models per lake
best_norm <- left_join(best_all, lake_meta,
                       by = c("lake" = "Lake.Short.Name")) |>
  left_join(kw) |> ungroup() |>
  mutate(across(c(4:24, 30:43),
                function(x)(x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))

filter(best_norm_a, best_met == "rmse") |>
  select(c(2, 28, 30:40, 43)) |> select(c(-1,-2,-9, -10)) |>
  cor() |> corrplot::corrplot()

disttance <- filter(best_norm_a, best_met == "rmse") |>
  select(c(2, 28, 30:40, 43)) |> select(c(-1,-2,-9, -10)) |>
  #mutate(model = factor(model, model, model),
  #       Reservoir.or.lake. = factor(Reservoir.or.lake.,
  #                                   Reservoir.or.lake.,
  #                                   Reservoir.or.lake.)) |>
  dist()
mydata.hclust <- hclust(disttance)

dat <- dendro_data(as.dendrogram(mydata.hclust), type = "rectangle")$segments
labs <- dendro_data(as.dendrogram(mydata.hclust), type = "rectangle")$labels |>
  arrange(as.numeric(label)) |> cbind(filter(best_norm_a, best_met == "rmse")[, 2:3])

ggplot() + geom_segment(data = dat,
                        aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data = labs, aes(x = x, y = y, label = lake, col = model),
            angle = -90, nudge_y = -1, hjust = 0) + theme_void() +
  theme(plot.margin = margin(b = 10)) + ylim(-15, 75)

# model as glm of the factors? 
  
##--------------- plots looking at the best performing parameter set -------------


# append lake meta data to the best parameter sets
best_rmse <- left_join(best_rmse, lake_meta, by = c("lake" = "Lake.Short.Name"))
best_r <- left_join(best_r, lake_meta, by = c("lake" = "Lake.Short.Name"))


# distribution of the best parameter set for each model
p_dist_lake_model <- best_rmse |> slice(which.min(rmse)) |>
  ggplot() + geom_histogram(aes(x = rmse, y = ..density.., fill = model),
                            position = 'stack', bins = 30, col = 1)  + 
  geom_density(aes(x = rmse, y =..density.., col = model), lwd = 1.5) +  
  thm + xlab("RMSE (K)") +
  facet_wrap(~model) + theme(legend.position = "none")


# distribution of the single best model per lake according to rmse
p_dist_lake <- best_rmse |> group_by(lake) |> filter(rmse == min(rmse)) |>
  ggplot() +
  geom_histogram(aes(x = rmse, y = ..density..),
                 bins = 20, col = 1) + 
  thm + xlab("RMSE (K)") +
  geom_density(aes(x = rmse, y =..density..), lwd = 1.5)

# pie chart with number each model performes best in a lake
p_pie <- best_rmse |> group_by(lake) |> slice(which.min(rmse)) |> group_by(model) |>
  summarise(n_best = n()) |> arrange(rev(model)) |> ggplot() +
  geom_bar(aes(x = "", y = n_best, fill = model), stat = "identity",
           col = 1) +
  geom_text(aes(x = "", y = n_best, label = n_best),
            position = position_stack(vjust=0.5), size = 5) +
  coord_polar("y", start = 0) + theme_void(base_size = 14) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(legend.position = "right")

# combine the two previous plots

png("best_fit_rmse.png", width = 11, height = 7, units = "in", res = 300)
subvp <- viewport(width = 0.35, height = 0.35, x = 0.8, y = 0.75)
p_dist_lake
print(p_pie, vp = subvp)
dev.off()


# distribution of the single best model per lake for all 6 metrics
p_dist_lake_a <- best_all_a |> pivot_longer(4:9) |> filter(best_met == name) |>
  ggplot() +
  geom_histogram(aes(x = value, y = ..count..),
                 bins = 20, col = 1) +
  thm + xlab("") +
  facet_wrap(~best_met, scales = "free")

ggsave("best_fit.pdf", p_dist_lake_a, width = 13, height = 7)

# a map of the lakes with the location color coded according to the best
# performing model
world <- map_data("world")
best_rmse |> group_by(lake) |> slice(which.min(rmse)) |> ggplot() +
  geom_map(
    data = world, map = world, fill = "grey33",
    aes(map_id = region)
  ) +
  geom_point(aes(y = latitude.dec.deg, x = longitude.dec.deg, col = model),
             size = 2, position = "jitter") + ylim(-40, 70) +
  xlim(-150, 180) + theme_void()

# a map of the lakes with the location color coded according to the lowest
# rmse from all four models
best_rmse |> group_by(lake) |> slice(which.min(rmse)) |> ggplot() +
  geom_map(
    data = world, map = world, fill = "grey33",
    aes(map_id = region)
  ) +
  geom_point(aes(y = latitude.dec.deg, x = longitude.dec.deg, col = rmse),
             size = 2, position = "jitter") + ylim(-90, 90) +
  xlim(-180, 180) + xlab("Longitude") + ylab("Latitude") +
  theme_pubclean(base_size = 17) + scale_color_viridis_c("RMSE (K)",option = "C")

ggsave("map_best_rmse.png", width = 13, height = 8)


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
  
  p_meta7 <- ggplot(data) + geom_point(aes_string(x = "Average.Secchi.disk.depth.m", y = measure, col = "model")) +
    ggtitle("min RMSE vs. Secchi disk depth") + thm + xlab("Average Secchi disk depth (m)") +
    ylab("RMSE (K)")
  
  p_meta5 <- ggplot(data) + geom_point(aes_string(x = "Duration", y = measure, col = "model")) +
    ggtitle("min RMSE vs.available data") + thm + xlab("Available observations (a)") +
    ylab("RMSE (K)")
  
  p_meta8 <- ggplot(data) + geom_point(aes_string(x = "months_meas", y = measure, col = "model")) +
    ggtitle("min RMSE vs. no. months with obs") + thm + xlab("Months with observations (-)") +
    ylab("RMSE (K)")
  
  p_meta9 <- ggplot(data) + geom_point(aes_string(x = "depth_meas", y = measure, col = "model")) +
    ggtitle("min RMSE vs. no unique depths") + thm + xlab("Unique depths with observations (-)") +
    ylab("RMSE (K)") + xlim(0, 200)
  
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
  
  p2_meta <- ggarrange(p_meta5,
                       p_meta8,
                       p_meta9,
                       p_meta10,
                       ncol = 2, nrow = 2, common.legend = TRUE)
  return(list(p1_meta, p2_meta))
}

# plot meta data vs best rmse for each lake and model
plot_meta(best_rmse)

# plot meta data vs best r for each lake and model
plot_meta(best_r, measure = "r")

# same plot but just for the best model per lake
best_rmse |> group_by(lake) |> filter(rmse == min(rmse)) |>  plot_meta()


##-------- compare models -----------------------

# distribution of rmse along all 2000 parameter sets and 73 lakes
res |> group_by(model) |> ggplot() +
  geom_violin(aes(x = model, y = rmse, fill = model)) + scale_y_log10() +
  thm

# check if the same models perform good or bad along the four models
best_rmse |> group_by(lake) |> reframe(range_rmse = range(rmse),
                                       sd_rmse = sd(rmse),
                                       min_rmse = min(rmse),
                                       max_rmse = max(rmse)) |>
  ggplot() + geom_col(aes(y = range_rmse, x = lake))

best_rmse |> group_by(lake) |> reframe(range_rmse = range(rmse),
                                       min_rmse = min(rmse),
                                       max_rmse = max(rmse),
                                       model_min = model[rmse == min(rmse)],
                                       model_max = model[rmse == max(rmse)]) |>
  ggplot() + geom_point(aes(x = min_rmse, y = max_rmse, col = model_max)) +
  geom_abline(aes(intercept = 0, slope = 1), col = 2, lty = 17)



##-------- relate calibrated parameter values to lake characteristics ----

# plot relating used wind speed scaling factor (of best rmse set) to lake area
# for each model
best_rmse |> ggplot() + geom_point(aes(y = wind_speed,
                                       x = lake.area.sqkm,
                                       col = model)) +
  thm + facet_wrap(.~model) + scale_x_log10()





##---------- plots for single models -----------------------------

# function to plot heatmaps for model performance along the different parameter
my_sens_plot <- function(m = "GLM", l = "Zurich", res_cali,
                         pars = c("swr", "wind_speed", "mixing.coef_mix_turb",
                                  "mixing.coef_mix_hyp"),
                         log = rep(FALSE, length(pars))) {
  
  spec <- viridis::viridis(15)
  dat <- filter(res_cali, model == m & lake == l)
  dat_best <- filter(dat, rmse < quantile(rmse, 0.1, na.rm = TRUE))
  thm <- theme_pubr(base_size = 13) + grids()
  
  pl <- lapply(combn(pars, 2, simplify = FALSE), function(p) {
    plt <- ggplot(dat_best) +
        geom_point(aes_string(x = p[1], y = p[2], color = "rmse"), shape = 15,
                   size = 2, alpha = 0) +
        geom_point(data = dat,
                   aes_string(x = p[1], y = p[2], color = "rmse"), shape = 15,
                   size = 1.5, alpha = 0.75) +
        scale_colour_gradientn(colours = rev(spec)) + thm +
        theme(legend.position = "none")
    
    if(log[pars %in% p[1]]) {
      plt <- plt + scale_x_log10()
    }
    if(log[pars %in% p[2]]) {
      plt <- plt + scale_y_log10()
    }
    return(plt)
  })
  
  t <- ggplot(dat) +
    geom_point(data = dat,
               aes(x = swr, y = wind_speed, color = rmse), shape = 15,
               size = 1.5, alpha = 0.75) +
    scale_colour_gradientn(colours = rev(spec))
  legend <- get_legend(t)
  
  t <- as_ggplot(arrangeGrob(text_grob(paste0("model: ", m, "\n lake: ", l)),
                             legend, ncol = 2))
  
  pl2 <- rep(list(NULL), (length(pars)-1)^2)
  
  k <- 1
  for (i in 1:(length(pars)-1)^2) {
    if(lower.tri(matrix(1:(length(pars)-1)^2,
                        ncol = (length(pars)-1),
                        nrow = (length(pars)-1),
                        byrow = TRUE),
                 diag = TRUE)[i]) {
      
      pl2[[i]] <- pl[[k]]
      k <- k+1
      
      if(!(i %in% c(1:(length(pars)-2),
                    length(pars)-1))) {
        pl2[[i]] <- pl2[[i]] + ylab("")
      }
      if(!(i %in% c(seq(length(pars)-1, (length(pars)-1)^2, by = length(pars)-1),
                    length(pars)-1))) {
        pl2[[i]] <- pl2[[i]] + xlab("")
      }
      if(i %in% seq(1, (length(pars)-1)^2, by = length(pars))) {
        pl2[[i]] <- ggMarginal(pl2[[i]], type = "densigram")
      }
      
    }
  }
  
  pl2[[(length(pars)-1)^2 - (length(pars)-1) + 1]] <- t
  
  do.call(grid.arrange, c(pl2, ncol = length(pars)-1, as.table = FALSE) )
  
  
  
}

# example plot showing a subset of parameter for lake Zurich and model GLM
my_sens_plot(res_cali = res)

# example plot showing all parameter for lake Zurich and model GLM
my_sens_plot(res_cali = res, pars = c("swr", "wind_speed", "Kw",
                                      "mixing.coef_mix_turb",
                                      "mixing.coef_mix_hyp",
                                      "mixing.coef_mix_conv"))

# example plot showing all parameter for lake Biel and model GOTM
my_sens_plot(res_cali = res, pars = c("swr", "wind_speed", "Kw",
                                      "turb_param.const_num",
                                      "bottom.h0b",
                                      "turb_param.k_min"),
             log = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
             l = "Biel", m = "GOTM")

# example plot showing all parameter for lake Kivu and model Simstrat
my_sens_plot(res_cali = res, pars = c("swr", "wind_speed", "Kw",
                                      "cd",
                                      "hgeo",
                                      "a_seiche"),
             l = "Kivu", m = "Simstrat")
