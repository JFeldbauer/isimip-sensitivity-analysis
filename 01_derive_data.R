## read in lake characteristics, calculate additional characteristics,
# do a cluster anqalysis and pick best perforing parameter sets

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

## if not exsiting create a folder for intermediate results
if(!file.exists("data_derived")) {
  dir.create("data_derived")
}

if(!file.exists("Plots")) {
  dir.create("Plots")
}
##-------------- read in data -------------------------------------------------

# source settings script
source("0_Settings.R")

# the results from the calibration runs descriptions of the different columns
# can be found in "data/results_lhc_description.csv"
res <- read.csv("data/results_lhc.csv")
# meta information about the lakes a description of the columns can be
# found in "data/Lake_meta_description.csv"
lake_meta <- read.csv("data/Lake_meta.csv")

lake_meta_desc <- read.csv("data/Lake_meta_description.csv")

## read in average lake water temperatures for the period 1980 -- 2000 from
# the climate simulations
lake_av_temp <- read.csv("data/av_temp_hist_sim.csv")
lake_av_temp_desc <- read.csv("data/av_temp_hist_sim_description.csv")

# read in hypsographic data for all lakes
# description of teh columns can be found in the file
# "data/lake_hypsographs_description.csv"
hyps <- read.csv("data/lake_hypsographs.csv")

##------------- pick best performing parameter sets -------------------------

# select best performing parameter set for each performance metric, lake,
# and model
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
  mutate(nmae = ifelse(nmae > 1e5, 0, nmae)) |>
  filter(nmae != 0 & !is.infinite(nmae)) |>
  slice_min(nmae) |> mutate(best_met = "nmae")

# extract best (mae) parameter set for each lake and model
best_mae <- res |> group_by(lake = lake,
                            model = model) |>
  slice_min(mae) |> mutate(best_met = "mae")

# data frame with all metrics for best set per lake, model, and metric
best_all <- rbind(best_bias, best_mae, best_nmae,
                  best_nse, best_r, best_rmse) |> ungroup()

# filter only the performance metrics which we set in 0_Settings.R
best_all <- filter(best_all, best_met %in% p_metrics) |>
  select(-which(colnames(best_all) %in% setdiff(c("bias",
                                                  "mae",
                                                  "nmae",
                                                  "nse",
                                                  "r",
                                                  "rmse"), p_metrics)))


saveRDS(best_all, "data_derived/best_par_sets.RDS")

## select single best performing model for each lake and performance metric
# extract best (rmse) parameter set for each lake
s_best_rmse <- res |> group_by(lake = lake) |>
  slice_min(rmse) |> mutate(best_met = "rmse")
# extract best (r) parameter set for each lake
s_best_r <- res |> group_by(lake = lake) |>
  slice_max(r) |> mutate(best_met = "r")

# extract best (nse) parameter set for each lake
s_best_nse <- res |> group_by(lake = lake) |>
  slice_max(nse) |> mutate(best_met = "nse")

# extract best (bias) parameter set for each lake
s_best_bias <- res |> group_by(lake = lake) |>
  slice_min(abs(bias)) |> mutate(best_met = "bias")

# extract best (nmae) parameter set for each lake
s_best_nmae <- res |> group_by(lake = lake) |>
  mutate(nmae = ifelse(nmae > 1e5, 0, nmae)) |>
  filter(nmae != 0 & !is.infinite(nmae)) |>
  slice_min(nmae, na_rm = TRUE) |> mutate(best_met = "nmae")

# extract best (mae) parameter set for each lake
s_best_mae <- res |> group_by(lake = lake) |>
  slice_min(mae) |> mutate(best_met = "mae")

# data frame with all metrics for single best model per lake
s_best_all <- rbind(s_best_bias, s_best_mae, s_best_nmae,
                    s_best_nse, s_best_r, s_best_rmse) |> ungroup()

# filter only the performance metrics which we set in 0_Settings.R
s_best_all <- filter(s_best_all, best_met %in% p_metrics) |>
  select(-which(colnames(s_best_all) %in% setdiff(c("bias",
                                                    "mae",
                                                    "nmae",
                                                    "nse",
                                                    "r",
                                                    "rmse"), p_metrics)))

saveRDS(s_best_all, "data_derived/single_best_model.RDS")

##---------------- gather additional lake characteristics --------------------

# get average Kw values from all best measures for each lake
kw <- best_all |> group_by(lake) |> reframe(kw = mean(Kw),
                                            kw_sd = sd(Kw)) |> ungroup()


# categorize lake by hypsography in three groups: convex, concave, neither
hyps_type <- hyps |> group_by(lake) |> mutate(area = area/max(area),
                                              level = 1 - depth/max(depth)) |>
  mutate(crv = case_when(mean(area - level < -0.025) >= 0.66 ~ "concave",
                         mean(area - level > 0.025) >= 0.66 ~ "convex",
                         .default = "neither")) |>
  mutate(crv = factor(crv, crv, crv)) |>
  select(lake, crv) |>
  group_by(lake) |> slice_head(n = 1) |> ungroup()


## calculate Osgood Index
lake_meta <- lake_meta |> mutate(osgood = mean.depth.m/(sqrt(mean.depth.m)))

##----------- cluster analysis ----------------------------------------------

# gather lake meta data for clustering
lake_meta <- lake_meta |> left_join(kw, by = c("Lake.Short.Name" = "lake")) |>
  left_join(hyps_type, by = c("Lake.Short.Name" = "lake")) |>
  left_join(lake_av_temp, by = c("Lake.Short.Name" = "lake")) |>
  mutate(Reservoir.or.lake. = factor(Reservoir.or.lake.)) 

# look at distribution of lake characteristics
lake_meta |> select(-Lake.Name, -Lake.Name.Folder, -Country,
                    -Average.Secchi.disk.depth.m,
                    -Light.extinction.coefficient.m, -kw_sd, -tsurf_sd,
                    - depth_meas, -tbot_sd) |>
  mutate(crv = as.numeric(crv),
         Reservoir.or.lake. = as.numeric(Reservoir.or.lake.)) |>
  pivot_longer(-1) |> ggplot() + geom_histogram(aes(x = value)) +
  facet_wrap(~name, scales = "free_x")

# correlation plot
lake_meta |> select(-Lake.Name, -Lake.Name.Folder, -Lake.Short.Name, -Country,
                    -crv, -Reservoir.or.lake., -Average.Secchi.disk.depth.m,
                    -Light.extinction.coefficient.m, -kw_sd, -tsurf_sd,
                    -tbot_sd) |>
  cor() |> corrplot::corrplot()
# pca with all data
lake_meta |> select(-Lake.Name, -Lake.Name.Folder, -Lake.Short.Name, -Country,
                    -Average.Secchi.disk.depth.m,
                    -Light.extinction.coefficient.m, -kw_sd, -tsurf_sd,
                    - depth_meas, -tbot_sd) |>
  mutate(across(contains(c("elevation.m",
                           "max.depth.m",
                           "mean.depth.m",
                           "lake.area.sqkm")),
                function(x)(log10(x - min(x) + 1)))) |> # log transform
  mutate(across(!contains(c("lake",
                           "Reservoir.or.lake.",
                           "crv")),
                function(x)(x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))) |>
  mutate(Reservoir.or.lake. = as.numeric(Reservoir.or.lake.),
         crv = as.numeric(crv)) |>
  prcomp() |> biplot()


# only select a subset of meta charcteristics
dat_clust <- lake_meta |>
  select(-Lake.Name, -Lake.Name.Folder, -Country, -Average.Secchi.disk.depth.m,
         -Light.extinction.coefficient.m, -months_median, -elevation.m,
         -depth_meas, -kw_sd, -tsurf_sd, -tbot_sd, -Reservoir.or.lake.,
         -tbot, -min_tsurf, -max.depth.m, -osgood, -crv) |>
  rename(lake = "Lake.Short.Name")
  

# sqrt transform  mean depth, max depth. lake area, and elevation and
# z-score normalize data
dat_clust_norm <- dat_clust |>
  mutate(across(contains(c("elevation.m",
                           "max.depth.m",
                           "mean.depth.m",
                           "lake.area.sqkm")),
                function(x)(log10(x - min(x) + 1)))) |> # log transfomr
  mutate(across(!contains(c("lake",
                            "Reservoir.or.lake.",
                            "crv")),
                function(x)(x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))


## kmeans clustering
# estimate optimal number of clusters using kmeans clustering
silhouette_score <- function(k, dat){
  km <- kmeans(dat, centers = k, nstart = 500, iter.max = 500)
  ss <- silhouette(km$cluster, dist(dat))
  mean(ss[, 3])
}
k <- 2:10
avg_sil <- sapply(k, silhouette_score, select(dat_clust_norm, -lake))
plot(k, type='b', avg_sil, xlab='Number of clusters',
     ylab='Average Silhouette Scores', frame=FALSE)
# number of clusters to create
nclust <- 5
# number of cluster for kmeans clustering
nkmclust <- nclust

km <- kmeans(select(dat_clust_norm, -lake), nkmclust, nstart = 500,
             iter.max = 500)

kmcluster <- km$cluster

## hirarchical clustering
# calculate distance for hirarchical clustering
disttance <- select(dat_clust_norm, -lake) |> dist()
hclus <- hclust(disttance)
hcluster <- cutree(hclus, k = nclust)

## attach cluster to the meta data frame and the data from clustering for plots
lake_meta <- lake_meta |> cbind(data.frame(kmcluster = factor(kmcluster),
                                           hcluster = factor(hcluster)))
dat_clust <- dat_clust |> cbind(data.frame(kmcluster = factor(kmcluster),
                                           hcluster = factor(hcluster)))

# add info about kw, osgood, and cluster to lake meta data descripton data.frame
lake_meta_desc <- mutate(lake_meta_desc,
                         column_name = ifelse(column_name == "Lake.Short.Name",
                                              "lake",
                                              column_name)) |>
  rbind(data.frame(column_name = c("kw", "crv", "osgood"),
                                   description = c("Average calibrated light extinction factor",
                                                   "Curvature type",
                                                   "Osgood index"),
                                   short_description = c("Kw",
                                                         "hyps.",
                                                         "Osgood index"),
                                   unit = c("m⁻¹", "-", "-"))) |>
  rbind(lake_av_temp_desc)

saveRDS(lake_meta, "data_derived/lake_meta_data_derived.RDS")
saveRDS(lake_meta_desc, "data_derived/lake_meta_desc_derived.RDS")
##------------ plots ----------------------------------------------------------

## hypsographs
# plot all hypsographs and their hypsographic type
p_hyps <- hyps |> left_join(hyps_type) |> group_by(lake) |>
  mutate(area = area/max(area),
         level = 1 - depth/max(depth)) |>ungroup() |>
  ggplot() + geom_line(aes(y = area, x = level, col = crv),
                       lwd = 1.5) +
  facet_wrap(~lake) + coord_flip() +
  geom_abline(aes(slope = 1, intercept = 0), col = "black", lty = 15) +
  thm


## dendogram for hirarchcal clustering
# extract data for plotting with ggplot2
dat <- dendro_data(as.dendrogram(hclus), type = "rectangle")$segments
labs <- dendro_data(as.dendrogram(hclus), type = "rectangle")$labels |>
  arrange(as.numeric(label)) |> cbind(select(dat_clust, lake))

# plot dendogram
p_tree <- ggplot() + geom_segment(data = dat,
                                  aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data = labs, aes(x = x, y = y, label = lake),
            angle = -90, nudge_y = -1, hjust = 0) + theme_void(base_size = 18) +
  theme(plot.margin = margin(b = 10)) + ylim(-6, 17) +
  geom_hline(aes(yintercept = 7.35), col = "grey42", lty = "dashed") +
  scale_color_viridis_d("best model", option = "H")



## PCA
pca_dat <- select(dat_clust_norm, -lake) |> prcomp()

p_pca <- as.data.frame(pca_dat$x) |>
  cbind(data.frame(kmcluster = factor(kmcluster),
                   hcluster = factor(hcluster))) |>
  ggplot() +
  geom_point(aes(x = PC1, y = PC2, col = kmcluster), size = 2.75) +
  geom_text(aes(x = PC1, y = PC2, col = kmcluster,
                label = dat_clust$lake),
            nudge_x = 0, nudge_y = 0.25) + 
  geom_segment(data = as.data.frame(pca_dat$rotation*10),
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm"))) + 
  geom_text(data = as.data.frame(pca_dat$rotation*10),
            aes(x = PC1, y = PC2, label = rownames(pca_dat$rotation)),
            col = "grey42") + thm +
  scale_color_viridis_d("Cluster") +
  xlab(paste0("PC1 ( ",
              round((pca_dat$sdev^2/sum(pca_dat$sdev^2))[1]*100, 1),
              "% )")) +
  ylab(paste0("PC1 ( ",
              round((pca_dat$sdev^2/sum(pca_dat$sdev^2))[2]*100, 1),
              "% )"))


p_bmc <- s_best_all |>
  left_join(dat_clust) |>
  select(lake, model, kmcluster, best_met) |>
  ggplot() + geom_histogram(aes(fill = model, x = kmcluster),
                            stat = "count", position = "Dodge") +
  thm +
  scale_fill_viridis_d("best model", option = "H") +
  facet_wrap(~best_met)

ggsave("Plots/best_model_per_clust.pdf", width = 13, height = 9)

p_rmsec <- s_best_all |>
  left_join(dat_clust) |>
  select(lake, model, kmcluster, !!p_metrics, best_met) |>
  pivot_longer(!!p_metrics) |> slice(which(best_met == name)) |>
  ggplot() + geom_violin(aes(y = value, x = kmcluster, fill = kmcluster)) +
  geom_jitter(aes(y = value,  x = kmcluster), height = 0,
              width = 0.125, size = 2.5, col = "grey42", alpha = 0.5) +
  thm + scale_fill_viridis_d("Cluster") +
  facet_wrap(~best_met, scales = "free_y") + theme(legend.position = "top") +
  xlab("Cluster") + ylab("")


ggsave("Plots/pca_cluster.pdf", p_pca, width = 11, height = 7, bg = "white")
ggsave("Plots/performance_cluster.pdf", p_rmsec, width = 11, height = 7, bg = "white")

# distributuin of the lake characteristics

p_clst_char <- lapply(colnames(dat_clust)[!colnames(dat_clust) %in% c("lake",
                                                                      "kmcluster",
                                                                      "hcluster")],
                      function(c) {
  dat <- select(dat_clust, c, "kmcluster")
  desc <- lake_meta_desc$short_description[lake_meta_desc$column_name == c]
  unit <- lake_meta_desc$unit[lake_meta_desc$column_name == c]
  if(is.factor(dat[, c])) {
    p <- dat |> table() |> as.data.frame() |> group_by(kmcluster) |>
       ggplot() +
      geom_col(aes_string(fill = c, y = "Freq", x = "kmcluster")) +
      scale_fill_viridis_d(desc, option = ifelse(c == "crv", "G", "E")) +
      xlab("") + thm +
      theme(legend.position = "top") + guides(fill=guide_legend(ncol=2))
  } else {
    p <- dat |> ggplot() +
      geom_violin(aes_string(y = c, fill = "kmcluster", x = "kmcluster")) +
      geom_jitter(aes_string(y = c,  x = "kmcluster"), height = 0,
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
  return(p)
  }) |> ggarrange(plotlist = _)


ggsave("Plots/clust_char.pdf", p_clst_char, width = 11, height = 7, bg = "white")


##### Cluster analysis and PCA, clustering by GOF -----

# z-score normalize data for single best model
best_norm_a <- left_join(s_best_all, lake_meta,
                         by = c("lake" = "Lake.Short.Name")) |>
  left_join(kw) |> left_join(hyps_type) |>
  mutate(osgood = mean.depth.m/(sqrt(mean.depth.m))) |>
  ungroup() |> mutate(across(c(32:35), function(x)(x - min(x) + 1e-4)^(1/3))) |>
  mutate(across(c(4:24, 30:43, 45),
                function(x)(x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))
# calculate distance for hirarchical clustering
distance2 <- filter(best_norm_a, best_met == "rmse") |>
  select(c(4, 6, 7, 8)) |>
  dist()
mydata.hclust2 <- hclust(distance2)

# extract data for plotting with ggplot2
dat2 <- dendro_data(as.dendrogram(mydata.hclust2), type = "rectangle")$segments
labs2 <- dendro_data(as.dendrogram(mydata.hclust2), type = "rectangle")$labels |>
  arrange(as.numeric(label)) |> cbind(filter(best_norm_a, best_met == "rmse")[, 2:3])

# plot dendogram
p_tree2 <- ggplot() + geom_segment(data = dat2,
                                   aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data = labs2, aes(x = x, y = y, label = lake, col = model),
            angle = -90, nudge_y = -1, hjust = 0) + theme_void(base_size = 18) +
  theme(plot.margin = margin(b = 10)) + ylim(-6, 12) +
  geom_hline(aes(yintercept = 5.5), col = "grey42", lty = "dashed") +
  scale_color_viridis_d("best model", option = "H")


clust2 <- cutree(mydata.hclust2, h = 5.5)

## PCA
pca_dat2 <- filter(best_norm_a, best_met == "rmse") |>
  select(c(4, 5, 6, 8, 10, 30:34, 40, 43:44, 46)) |>
  mutate(crv = as.numeric(crv)) |> prcomp()

plot(pca_dat2)
biplot(pca_dat2)

p_pca2 <- as.data.frame(pca_dat2$x) |> cbind(clust2) |>
  mutate(clust = factor(clust2)) |> ggplot() +
  geom_point(aes(x = PC1, y = PC2, col = clust), size = 2.75) +
  geom_text(aes(x = PC1, y = PC2, col = clust,
                label = dat_clust$lake),
            nudge_x = 0, nudge_y = 0.5) + 
  geom_segment(data = as.data.frame(pca_dat2$rotation*10),
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm"))) + 
  geom_text(data = as.data.frame(pca_dat2$rotation*10),
            aes(x = PC1, y = PC2, label = rownames(pca_dat2$rotation)),
            col = "grey42") + theme_minimal(base_size = 18) +
  scale_color_viridis_d("Cluster") +
  xlab(paste0("PC1 ( ",
              round((pca_dat2$sdev^2/sum(pca_dat2$sdev^2))[1]*100, 1),
              "% )")) +
  ylab(paste0("PC1 ( ",
              round((pca_dat2$sdev^2/sum(pca_dat2$sdev^2))[2]*100, 1),
              "% )"))
