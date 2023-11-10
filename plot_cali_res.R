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
  filter(nmae != 0 & !is.infinite(nmae)) |>
  slice_min(nmae, na_rm = TRUE) |> mutate(best_met = "nmae")

# extract best (mae) parameter set for each lake
best_mae_a <- res |> group_by(lake = lake) |>
  filter(nmae != 0 & !is.infinite(nmae)) |>
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

# look at hypsographs and try to use them to classify lakes
hyps <- read.csv("data/lake_hypsographs.csv")

# categorize lake by hypsography in three groups: convex, concave, neither
hyps_type <- hyps |> group_by(lake) |> mutate(area = area/max(area),
                                              level = 1 - depth/max(depth)) |>
  mutate(crv = case_when(mean(area - level < -0.025) >= 0.66 ~ "concave",
                         mean(area - level > 0.025) >= 0.66 ~ "convex",
                         .default = "neither")) |>
  mutate(crv = factor(crv, crv, crv)) |>
  select(lake, crv) |>
  group_by(lake) |> slice_head(n = 1) |> ungroup()

# plot all hypsographs and their hypsographic type
hyps |> left_join(hyps_type) |> group_by(lake) |>
  mutate(area = area/max(area),
         level = 1 - depth/max(depth)) |>ungroup() |>
  ggplot() + geom_line(aes(y = area, x = level, col = crv),
                       lwd = 1.5) +
  facet_wrap(~lake) + coord_flip() +
  geom_abline(aes(slope = 1, intercept = 0), col = "black", lty = 15) +
  theme_minimal(base_size = 14)

ggsave("hypso_class.png", width = 13, height = 13, bg = "white")

# z-score normalize data for single best model
best_norm_a <- left_join(best_all_a, lake_meta,
                         by = c("lake" = "Lake.Short.Name")) |>
  left_join(kw) |> left_join(hyps_type) |>
  ungroup() |> mutate(across(c(32:35), function(x)(x - min(x) + 1e-4)^(1/3))) |>
  mutate(across(c(4:24, 30:43),
                function(x)(x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))

# calculate distance for hirarchical clustering
disttance <- filter(best_norm_a, best_met == "rmse") |>
  select(c(2, 28, 30:40, 43:44)) |> select(c(-1,-2,-9, -10, -11)) |>
  dist()
mydata.hclust <- hclust(disttance)

# extract data for plotting with ggplot2
dat <- dendro_data(as.dendrogram(mydata.hclust), type = "rectangle")$segments
labs <- dendro_data(as.dendrogram(mydata.hclust), type = "rectangle")$labels |>
  arrange(as.numeric(label)) |> cbind(filter(best_norm_a, best_met == "rmse")[, 2:3])

# plot dendogram
p_tree <- ggplot() + geom_segment(data = dat,
                        aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data = labs, aes(x = x, y = y, label = lake, col = model),
            angle = -90, nudge_y = -1, hjust = 0) + theme_void(base_size = 18) +
  theme(plot.margin = margin(b = 10)) + ylim(-6, 12) +
  geom_hline(aes(yintercept = 5.5), col = "grey42", lty = "dashed") +
  scale_color_viridis_d("best model", option = "H")


clust <- cutree(mydata.hclust, h = 5.5)

## PCA
pca_dat <- filter(best_norm_a, best_met == "rmse") |>
  select(c(2, 28, 30:40, 43:44)) |> select(c(-1,-2,-9, -10, -11)) |>
  mutate(crv = as.numeric(crv)) |> prcomp()

plot(pca_dat)
biplot(pca_dat)

p_pca <- as.data.frame(pca_dat$x) |> cbind(clust) |>
  mutate(clust = factor(clust)) |> ggplot() +
  geom_point(aes(x = PC1, y = PC2, col = clust), size = 2.75) +
  geom_text(aes(x = PC1, y = PC2, col = clust,
                label = best_bias_a$lake),
            nudge_x = 0, nudge_y = 1.5) + 
  geom_segment(data = as.data.frame(pca_dat$rotation*10),
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm"))) + 
  geom_text(data = as.data.frame(pca_dat$rotation*10),
            aes(x = PC1, y = PC2, label = rownames(pca_dat$rotation)),
            col = "grey42") + theme_minimal(base_size = 18) +
  scale_color_viridis_d("Cluster") +
  xlab(paste0("PC1 ( ",
              round((pca_dat$sdev^2/sum(pca_dat$sdev^2))[1]*100, 1),
              "% )")) +
  ylab(paste0("PC1 ( ",
              round((pca_dat$sdev^2/sum(pca_dat$sdev^2))[2]*100, 1),
              "% )"))

p_bmc <- best_all_a |>
  left_join(data.frame(lake = best_rmse_a$lake,
                       cluster = factor(clust))) |>
  select(lake, model, cluster, best_met) |>
  ggplot() + geom_histogram(aes(fill = model, x = cluster),
                            stat = "count", position = "Dodge") +
  theme_minimal(base_size = 18) +
  scale_fill_viridis_d("best model", option = "H") +
  facet_wrap(~best_met)

p_rmsec <- best_all_a |>
  left_join(data.frame(lake = best_rmse_a$lake,
                       cluster = factor(clust))) |>
  select(lake, model, cluster, best_met, rmse, nse, r, bias, mae, nmae) |>
  mutate(nmae = ifelse(nmae > 1e4, NA, nmae)) |>
  pivot_longer(5:10) |> slice(which(best_met == name)) |>
  ggplot() + geom_boxplot(aes(y = value, x = cluster, fill = cluster)) +
  theme_minimal(base_size = 18) + scale_fill_viridis_d("Cluster") +
  facet_wrap(~best_met, scales = "free_y")


ggarrange(p_tree, p_pca, p_bmc, p_rmsec, nrow = 2, ncol = 2, labels = "AUTO")

ggsave("clustering_sbest.png", width = 20, height = 20, bg = "white")

# distributuin of the lake characteristics

lake_meta |> select(-1, -(3:5), -12, -13, -14, -17, -18) |>
  rename(lake = "Lake.Short.Name") |> left_join(kw) |>
  left_join(select(mutate(hyps_type, curvature = as.numeric(crv)), -crv)) |>
  left_join(data.frame(lake = best_rmse_a$lake,
                       cluster = factor(clust, clust, clust))) |>
  pivot_longer(2:11) |>
  ggplot() +
  geom_boxplot(aes(x = cluster, y = value, fill = cluster)) +
  scale_fill_viridis_d("Cluster") +
  facet_wrap(~name, scales = "free") + theme_minimal(base_size = 18)

ggsave("clust_char.png", width = 13, height = 9, bg = "white")

##--------------- statistical models for best model/performacne ----------------
# fit multinomial log-linear model 
dat <- left_join(best_rmse_a, lake_meta,
                 by = c("lake" = "Lake.Short.Name")) |>
  left_join(kw) |> left_join(hyps_type) |>
  mutate(model = factor(model, model, model),
         Reservoir.or.lake. = factor(Reservoir.or.lake.,
         Reservoir.or.lake.,
         Reservoir.or.lake.)) |> ungroup()

dat$model <- relevel(dat$model, ref = "Simstrat")
test <- multinom(model ~ kw + elevation.m + max.depth.m +
                   lake.area.sqkm + latitude.dec.deg + longitude.dec.deg + crv,
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
                    crv = levels(dat$crv))
  

res_m <- predict(test, newdata = pdat, "probs")
res_m <- cbind(res_m, pdat) |> pivot_longer(1:4)

res_m |> ggplot() + geom_line(aes(x = lake.area.sqkm, y = value, col = name)) +
  facet_grid(crv~kw)

## estimate performance metric based upon lake characteristics

ggplot(dat) + geom_point(aes(x = kw, y = lake.area.sqkm, col = model, pch = crv)) +
  scale_y_log10() + scale_x_log10()


## model best rsme from lake propertiers

rmse_m <- left_join(best_all_a, lake_meta,
                    by = c("lake" = "Lake.Short.Name")) |>
  left_join(kw) |> left_join(hyps_type) |> lm(formula = rmse ~ (kw + elevation.m + max.depth.m +
               lake.area.sqkm + latitude.dec.deg + longitude.dec.deg +
               reldepth_median + months_meas) * crv)

summary(rmse_m)
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
