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
library(akima)

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

##----------------- first plots ---------------------------------

# plot delta sensitivity metrics for one lake
res |> filter(lake == "Annie" & var == "rmse") |> ggplot() +
  geom_col(aes(y = delta, x = names)) +
  geom_errorbar(aes(x = names, ymin = delta - delta_conf/2,
                    ymax = delta + delta_conf), col = 2) +
  theme_pubr() + grids() + facet_wrap(~model, scales = "free_x") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  xlab("parameter") + ggtitle("Annie")


# correlation between S1 and delta
res_o |> ggplot(aes(x = delta, y = S1, col = model), size = 0.6,
                alpha = 0.666) +
  geom_abline(aes(intercept = 0, slope = 1), col = "grey42", lty = "dashed") +
  geom_point() +
  facet_grid(var ~ names) + theme_pubr(base_size = 17) + grids() +
  xlim(0, 1) +  ylim(0, 1)

##-------------- interaction measure ---------------------------------------

# interactions measure: S_interact = 1 - sum(S1_i)
dat_iat <- res_o |> group_by(lake, model, var) |> reframe(iat = 1 - sum(S1)) |>
  left_join(meta, by = c("lake" = "Lake.Short.Name"))

dat_iat |> 
  mutate(var = ifelse(var == "bias", var, toupper(var))) |>
  ggplot() + geom_boxplot(aes(x = model,
                              y = iat, fill = var)) +
  thm + xlab("Model") + ylab("Interaction measure") +
  scale_fill_viridis_d("Performance metric", option = "H")

dat_iat <- res_o |> group_by(lake, model, var) |> reframe(iat = 1 - sum(S1)) |>
  left_join(meta, by = c("lake" = "Lake.Short.Name"))
# same plot but split to clusters
dat_iat |> 
  mutate(var = ifelse(var == "bias", var, toupper(var))) |>
  ggplot() + geom_boxplot(aes(x = model,
                              y = iat, fill = var)) +
  facet_wrap(~kmcluster) +
  thm + xlab("Model") + ylab("Interaction measure") +
  scale_fill_viridis_d("Performance metric", option = "H")

ggsave("Plots/interaction_clust.pdf", width = 11, height = 8)

dat_iat |> filter(iat > 0.2) |> ggplot() +
  geom_histogram(aes(x = 1, fill = model, y = after_stat(count/sum(count)))) +
  scale_fill_viridis_d("Model", option = "C", end = 0.9) + thm

dat_iat |> filter(iat > 0.2) |> ggplot() +
  geom_histogram(aes(x = var, fill = model), stat = "count") +
  scale_fill_viridis_d("Model", option = "C", end = 0.9) +
  facet_wrap(.~kmcluster) + thm

# compare iat with observed data
dat_iat |>
  ggplot() + geom_point(aes(x = Duration, y = iat, col = kmcluster), size = 3) +
  facet_grid(~model) + xlab("Years with observations (-)") +
  ylab("Interaction measure (-)") +
  thm + scale_color_viridis_d("Cluster")

ggsave("Plots/iat_obs_duration.pdf", width = 13, height = 5)

# relate iat to most sensitive parameter
res_o |> group_by(lake, model, var) |> reframe(max_d = max(delta),
                                               max_S1 = max(S1)) |>
  left_join(dat_iat) |> ggplot() +
  geom_point(aes(iat, max_d, col = kmcluster)) + facet_grid(.~model) +
  scale_color_viridis_d("Cluster") + thm

##--------- single most sensitive parameter  ---------------------------

# boxplots of sensitivity metrics
res_o |> pivot_longer(cols = c(delta, S1)) |> filter(names != "dummy") |>
  mutate(name = case_match(name,
                           "delta" ~ "\u03B4",
                           "S1" ~ "S1")) |>
  mutate(var = ifelse(var == "bias", var, toupper(var))) |>
  ggplot() +
  geom_boxplot(aes(x = names, y = value, fill = name)) +
  facet_grid(var ~ model, scales = "free") +
  scale_fill_manual("Sensitivity \n measure",
                    values = c("#45B2DD", "#72035F")) +
  thm + xlab("parameter") +
  theme(axis.text.x=element_text(angle = -60, hjust = 0)) 

ggsave("Plots/sens_value.pdf", width = 11, height = 9, device = cairo_pdf)


# calculate single most sensitive parameter for each lake and model
res_mip <- res |> group_by(lake, model, var) |>
  reframe(par_d = names[delta == max(delta)],
          par_S1 = names[S1 == max(S1)]) |> select(-par_S1)

# plot distribution of single most sensitive parameter over all lakes and models
# for all measures
res_mip |> pivot_longer(cols = par_d) |>
  mutate(name = case_match(name,
                           "par_d" ~ "\u03B4",
                           "par_S1" ~ "S1")) |> ggplot() +
  geom_histogram(aes(x = value, fill = name), stat = "count",
                 position = "dodge", col = 1) +
  facet_grid(var~model, scales = "free_x") + thm +
  grids() + xlab("parameter") +
  theme(axis.text.x=element_text(angle = -60, hjust = 0),
        legend.position = "none") +
  scale_fill_manual("", values = c("#72035F"))

ggsave("Plots/sens_single.png", width = 11, height = 9)


# look at cluster
res_mip |> pivot_longer(cols = par_d) |>
  mutate(name = case_match(name,
                           "par_d" ~ "\u03B4",
                           "par_S1" ~ "S1")) |>
  left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(var = ifelse(var == "bias", var, toupper(var))) |>
  ggplot() +
  geom_histogram(aes(x = value, fill = kmcluster),
                 stat = "count", position = "dodge", col = 1) +
  facet_grid(var~model, scales = "free_x") + thm +
  grids() + xlab("parameter") +
  theme(axis.text.x=element_text(angle = -60, hjust = 0)) +
  scale_fill_viridis_d("Cluster")


ggsave("Plots/sens_single_clust.png", width = 11, height = 9)


##--------- group of most sensitive parameters -------------------------------

# calculate group of most sensitive parameters for each lake and model
# function to return all values in a vector that contribute to frac (default 75)
# percent of the sum of all values
gmiv <- function(x, frac = .75) {
  tmp <- sort(x, decreasing = TRUE)

  
  id <- which(cumsum(tmp)/sum(tmp) < frac)
  if(length(id) > 0) {
    im <- (max(id)) + 1
  } else {
    im <- 1
  }

  id <- 1:im
  return(tmp[id])
}

# for each model, lake, and measure return the most sensitive parameters
res_gip <- res |> group_by(lake, model, var) |> 
  reframe(par_d_05 = paste0(names[delta %in% gmiv(delta, .5)], collapse = ", "),
          n_par_d_05 = length(names[delta %in% gmiv(delta, .5)]),
          par_S1_05 = paste0(names[S1 %in% gmiv(S1, .5)], collapse = ", "),
          n_par_S1_05 = length(names[S1 %in% gmiv(S1, .5)]),
          par_d_75 = paste0(names[delta %in% gmiv(delta, .75)], collapse = ", "),
          n_par_d_75 = length(names[delta %in% gmiv(delta, .75)]),
          par_S1_75 = paste0(names[S1 %in% gmiv(S1, .75)], collapse = ", "),
          n_par_S1_75 = length(names[S1 %in% gmiv(S1, .75)]),
          par_d_1 = paste0(names[delta %in% gmiv(delta, 1)], collapse = ", "),
          n_par_d_1 = length(names[delta %in% gmiv(delta, 1)]),
          par_S1_1 = paste0(names[S1 %in% gmiv(S1, 1)], collapse = ", "),
          n_par_S1_1 = length(names[S1 %in% gmiv(S1, 1)]))

# plot distribution of most sensitive parameters over all lakes and models
# for both measures
delta_gip <- res_gip |> pivot_longer(seq(4, 12, by = 4)) |>
  filter(grepl("par_d_.*", name)) |>
  mutate(frac = gsub(".*_", "", name)) |> 
  group_by(model, var, frac, lake) |> reframe(par = unlist(strsplit(value, ", ")),
                                                   meas = "\u03B4")
S1_gip <- res_gip |> pivot_longer(seq(6, 14, by = 4)) |>
   filter(grepl("par_S1_.*", name)) |>
   mutate(frac = gsub(".*_", "", name)) |> 
   group_by(model, var, frac, lake) |> reframe(par = unlist(strsplit(value, ", ")),
                                               meas = "S1")

#rbind(delta_gip, S1_gip) 
delta_gip |> filter(frac == "1") |>
  mutate(var = ifelse(var == "bias", var, toupper(var))) |>
  ggplot() +
  geom_histogram(aes(x = par, fill = meas), stat = "count", position = "dodge") +
  facet_grid(var~model, scales = "free_x") + thm +
  grids() + xlab("parameter") +
  theme(axis.text.x=element_text(angle = -55, hjust = 0),
        legend.position = "none") +
  scale_fill_manual("Sensitivity \n measure", values = c("#72035F"))

ggsave("Plots/count_sens.png", width = 14, height = 11)

# also plot for the different cluster
#rbind(delta_gip, S1_gip) |>
 delta_gip |> filter(frac == "1") |>
  left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(var = ifelse(var == "bias", var, toupper(var))) |>
  ggplot() + geom_histogram(aes(x = par,  fill = kmcluster),
                            stat = "count", position = "dodge",
                            col = 1) + 
  facet_grid(var~model, scales = "free_x") + thm + grids() +
  scale_fill_viridis_d("Cluster") +
  theme(axis.text.x=element_text(angle = -55, hjust = 0)) +
  xlab("Parameter")

ggsave("Plots/count_sens_clust.pdf", width = 14, height = 11)

# different plot with frequencies
#rbind(delta_gip, S1_gip)
delta_gip |> filter(frac == "1") |>
  left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  group_by(model, kmcluster, var, par) |> reframe(cnt = length(var)) |>
  group_by(model, kmcluster, var) |> mutate(cnt = cnt/sum(cnt)) |>
  ungroup() |>
  #complete(model, kmcluster, var, par) |>
  mutate(var = ifelse(var == "bias", var, toupper(var))) |>
  ggplot() + geom_tile(aes(x = par,  y = kmcluster, fill = cnt)) + 
  facet_grid(var~model, scales = "free_x") + thm + grids() +
  scale_fill_viridis_c("Frequency", option = "C") +
  theme(axis.text.x=element_text(angle = -55, hjust = 0)) +
  xlab("Parameter") + ylab("Cluster")

ggsave("Plots/freq_sens_clust.pdf", width = 14, height = 11)


## plot distribution of number of sensitive parameters
# rbind(delta_gip, S1_gip)
delta_gip |> group_by(model, var, frac, lake, meas) |>
  reframe(n = n()) |> mutate(frac = factor(frac, levels = c("1", "75", "05"),
                                           labels = c("100%", "75%", "50%"))) |>
  mutate(var = ifelse(var == "bias", var, toupper(var))) |>
  ggplot() + geom_histogram(aes(x = n,  fill = frac),
                            stat = "count", position = "dodge") + 
  facet_grid(var~model) + thm + grids() +
  scale_fill_viridis_d("Fraction of total sum", option = "D") +
  xlab("Number of parameters contributing")

ggsave("Plots/count_imp_par.pdf", width = 14, height = 11)

## alternative plot
# rbind(delta_gip, S1_gip)
delta_gip |> group_by(model, var, frac, lake, meas) |>
  reframe(n = n()) |> mutate(frac = factor(frac, levels = c("1", "75", "05"),
                                           labels = c("100%", "75%", "50%"))) |>
  group_by(model, frac, var, n) |> reframe(cnt = length(var)) |>
  complete(model, frac, var, n) |>
  ggplot() + geom_tile(aes(x = n,  y = frac, fill = cnt)) + 
  facet_grid(var~model) + theme_pubr(base_size = 17) + grids() +
  scale_fill_viridis_c("count", option = "C") +
  ylab("Fraction of total sum") +
  xlab("Number of parameters contributing")

ggsave("Plots/count_imp_par2.png", width = 14, height = 11)

## relate number of sensitive parameters to interaction measure
rbind(delta_gip, S1_gip) |> filter(frac == "1") |> group_by(model, var, frac, lake, meas) |> 
  reframe(n = n()) |> #pivot_wider(names_from = meas, values_from = n) |>
  left_join(dat_iat) |> #mutate(d_n = δ - S1) |> 
  ggplot() + geom_point(aes(x = n, y = iat, col = kmcluster), size = 3) + 
  facet_grid(var~model) + scale_color_viridis_d("Cluster") + thm


# lok at number of sensitive parameters, maximum sensitivity value and iat
res_o |> group_by(lake, model, var) |> reframe(max_d = max(delta),
                                               max_S1 = max(S1)) |>
  left_join(dat_iat) |> left_join(delta_gip |> filter(frac == "1") |> group_by(model, var, lake) |> 
                                    reframe(n = n())) |>
  ggplot() + geom_point(aes(x = iat, y = max_d, col = kmcluster), size = 2.5) +
  facet_grid(n ~ model) + scale_color_viridis_d("Cluster") + thm

# plot iat against metric value of best performing parameter set
rbind(S1_gip) |> filter(frac == "1") |> group_by(model, var, frac, lake, meas) |> 
  reframe(np = n()) |> left_join(best_all, by = c("lake" = "lake",
                                                 "model" = "model",
                                                 "var" = "best_met")) |>
  left_join(dat_iat, by = c("lake" = "lake",
                            "model" = "model",
                            "var" = "var")) |>
  select(lake, model, var, iat, np, kmcluster, rmse, r, bias, nse) |>
  pivot_longer(cols = 7:10) |> filter(var == name) |> select(-name) |>
  ggplot() + geom_point(aes(y = iat, x = value)) +
  facet_grid(.~var, scales = "free_x") +
  thm

##--------------- look at sensitivity of scaling factors -----------


# look at different cluster wind speed
res |> left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(delta = ifelse(delta > 1, 1, delta)) |>
  filter(names == "wind_speed") |> ggplot() +
  geom_boxplot(aes(x = as.numeric(kmcluster), y = delta, fill = kmcluster)) +
  facet_grid(var~model) + thm +
  scale_fill_viridis_d("Cluster") + xlab("Cluster") +
  ylab("Delta sensitivity wind scaling")

ggsave("Plots/sensitivity_wind_clust_model.png", width = 13, height = 9)

res |> left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(delta = ifelse(delta > 1, 1, delta)) |>
  filter(names == "wind_speed") |> ggplot() +
  geom_boxplot(aes(x = model, y = delta, fill = model)) +
  facet_grid(var~kmcluster) + thm +
  scale_fill_viridis_d("Model", option = "C", end = 0.9) + xlab("Cluster") +
  ylab("Delta sensitivity wind scaling") +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

ggsave("Plots/sensitivity_wind_clust.png", width = 13, height = 9)

# look at different cluster swr
res |> left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(delta = ifelse(delta > 1, 1, delta)) |>
  filter(names == "swr") |> ggplot() +
  geom_boxplot(aes(x = as.numeric(kmcluster), y = delta, fill = kmcluster)) +
  facet_grid(var~model) + thm +
  scale_fill_viridis_d("Cluster") + xlab("Cluster") +
  ylab("Delta sensitivity swr scaling")

ggsave("Plots/sensitivity_swr_clust_model.png", width = 13, height = 9)

res |> left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(delta = ifelse(delta > 1, 1, delta)) |>
  filter(names == "swr") |> ggplot() +
  geom_boxplot(aes(x = model, y = delta, fill = model)) +
  facet_grid(var~kmcluster) + thm +
  scale_fill_viridis_d("Model", option = "C", end = 0.9) + xlab("Cluster") +
  ylab("Delta sensitivity swr scaling") +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

ggsave("Plots/sensitivity_swr_clust.png", width = 13, height = 9)

# look at different cluster Kw (use same plot as for calibrated Kw values
# per cluster in script 02_)
res |> left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(delta = ifelse(delta > 1, 1, delta)) |>
  filter(names == "Kw") |> ggplot() +
  geom_boxplot(aes(x = as.numeric(kmcluster), y = delta, fill = kmcluster)) +
  facet_grid(var~model) + thm +
  scale_fill_viridis_d("Cluster") + xlab("Cluster") +
  ylab("Delta sensitivity Kw scaling")

ggsave("Plots/sensitivity_Kw_clust_model.png", width = 13, height = 9)

res |> left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(delta = ifelse(delta > 1, 1, delta)) |>
  filter(names == "Kw") |> ggplot() +
  geom_boxplot(aes(x = model, y = delta, fill = model)) +
  facet_grid(var~kmcluster) + thm +
  scale_fill_viridis_d("Model", option = "C", end = 0.9) + xlab("Cluster") +
  ylab("Delta sensitivity Kw scaling") +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

ggsave("Plots/sensitivity_Kw_clust.png", width = 13, height = 9)

##---------- plots for single models -----------------------------

#visual comparison of heatmaps against the calculated sensitivity metrics

# function to plot heatmaps for model performance along the different parameter
my_sens_plot <- function(m = "GLM", l = "Zurich", res_cali, res_sens,
                         smet = "rmse", contour = FALSE, n_contour = 5) {
  
  sens <- res_sens |> filter(lake == l & model == m & var == smet)
  pars <- unique(sens$names)
  pars_s <- pars
  pars <- pars[pars!="dummy"]
  log <- rep(FALSE, length(pars))
  log[pars == "turb_param.k_min"] <- TRUE
  
  

  dat <- filter(res_cali, model == m & lake == l)
  # calculate interaction quantity
  S_interaction <-  1 - sum(sens$S1)
  if (smet == "rmse") {
    dat_best <- filter(dat, rmse < quantile(rmse, 0.1, na.rm = TRUE))
  }
  if (smet == "r") {
    dat_best <- filter(dat, r > quantile(r, 0.9, na.rm = TRUE))
  }
  if (smet == "nse") {
    dat_best <- filter(dat, nse > quantile(nse, 0.9, na.rm = TRUE))
  }
  if (smet == "bias") {
    dat_best <- filter(dat, abs(bias) < quantile(abs(bias), 0.1, na.rm = TRUE))
  }
  thm <- theme_pubr(base_size = 13) + grids()
  
  pl <- lapply(combn(pars, 2, simplify = FALSE), function(p) {
    # interpolate for contour plot
    dat_tmp <- dat |> select(all_of(p[1]), all_of(p[2]), all_of(smet)) |>
      setNames(c("x", "y", "z"))
    xm <- mean(dat_tmp$x)
    ym <- mean(dat_tmp$y)
    dat_tmp <- mutate(dat_tmp, x = x/xm, y = y/ym)
    dat_tmp <- with(dat_tmp, interp(x = x, y = y, z = z, nx = 100, ny = 100, extrap = FALSE))
    cnt <- data.frame(dat_tmp$z)
    colnames(cnt) <- dat_tmp$y
    cnt$x <- dat_tmp$x

    cnt <- pivot_longer(cnt, cols = -x, names_to = "y") |>
      mutate(y = as.numeric(y)) |> mutate(x = x*xm, y = y*ym)
    
    # plot
    plt <- ggplot(dat_best) +
      geom_point(aes_string(x = p[1], y = p[2], color = smet), shape = 15,
                 size = 2, alpha = 0) +
      geom_point(data = dat,
                 aes_string(x = p[1], y = p[2], color = smet), shape = 15,
                 size = 1.5, alpha = 0.75) +
      scale_colour_viridis_c() + thm +
      theme(legend.position = "none",
            plot.margin = margin(t = 20,    # Top margin
                                 r = 30,    # Right margin
                                 b = 10,    # Bottom margin
                                 l = 10))
      
    if(contour) {
      plt <- plt + geom_contour(data = cnt, aes(x = x, y = y, z = value),
                                col = alpha("black", 0.5), lwd = 0.85, bins = n_contour)
    }
    
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
               aes_string(x = "swr", y = "wind_speed", color = smet), shape = 15,
               size = 1.5, alpha = 0.75) +
    scale_colour_viridis_c()
  legend <- get_legend(t)
  
  t <- as_ggplot(arrangeGrob(text_grob(paste0("model: ", m,
                                              "\n lake: ", l,
                                              "\n metric: ", smet)),
                             legend, ncol = 2))
  
  pl2 <- rep(list(NULL), (length(pars)-1)^2)
  
  ps1 <- sens |> filter(names != "dummy") |> ggplot() +
    geom_col(aes(y = delta, x = names), fill = "#45B2DD") +
    geom_errorbar(aes(x = names, ymin = delta - delta_conf/2,
                      ymax = delta + delta_conf),
                  col = "#0D3D20", lwd = 1.25, width=.5) +
    geom_hline(data = filter(sens, names == "dummy"),
               aes(yintercept = delta), col = "grey42",
               lty = "dashed", linewidth = 1.25) + 
    geom_hline(data = filter(sens, names == "dummy"),
               aes(yintercept = delta + delta_conf), col = "grey42",
               lty = "dashed", linewidth = 1) + 
    geom_hline(data = filter(sens, names == "dummy"),
               aes(yintercept = delta - delta_conf), col = "grey42",
               lty = "dashed", linewidth = 1) + 
    theme_pubr() + grids() +
    theme(axis.text.x=element_text(angle = -65,
                                   hjust = 0, size = 9),
          plot.margin = margin(t = -40,    # Top margin
                               r = 10,    # Right margin
                               b = -70,    # Bottom margin
                               l = 10)) + # Left margin
    xlab("") + ylab("")
    
  ps2 <- sens |> filter(names != "dummy") |> ggplot() +
      geom_col(aes(y = S1, x = names), fill = "#72035F") +
      geom_errorbar(aes(x = names, ymin = S1 - S1_conf/2,
                        ymax = S1 + S1_conf),
                    col = "#03DE67", lwd = 1.25, width=.5) +
    geom_hline(data = filter(sens, names == "dummy"),
               aes(yintercept = S1), col = "grey42",
               lty = "dashed", linewidth = 1.25) + 
    geom_hline(data = filter(sens, names == "dummy"),
               aes(yintercept = S1 + S1_conf), col = "grey42",
               lty = "dashed", linewidth = 1) + 
    geom_hline(data = filter(sens, names == "dummy"),
               aes(yintercept = S1 - S1_conf), col = "grey42",
               lty = "dashed", linewidth = 1) + 
      theme_pubr() + grids() +
    theme(axis.text.x=element_text(angle = -65,
                                   hjust = 0, size = 9),
          plot.margin = margin(t = -40,    # Top margin
                               r = 10,    # Right margin
                               b = -70,    # Bottom margin
                               l = 10)) + # Left margin
      xlab("") + ylab("")
  
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


  
  pl2[[21]] <- as_ggplot(text_grob("\u03B4"))
  pl2[[16]] <- as_ggplot(text_grob("S1"))
  pl2[[22]] <- ps1
  pl2[[17]] <- ps2
  pl2[[6]] <- t
  pl2[[11]] <- as_ggplot(text_grob(paste0("S_interaction = ",
                                          round(S_interaction, 3))))
  
  do.call(grid.arrange, c(pl2, ncol = length(pars)-1, as.table = FALSE) )
  
}

my_sens_plot(m = "GOTM", l = "Kivu", res_cali = res_cali,
                           res_sens = res_o)
ggsave("Plots/GOTM_kivu.png", width = 17, height = 12)

my_sens_plot(m = "GLM", l = "Biel", res_cali = res_cali,
                           res_sens = res_o)
ggsave("Plots/GLM_biel.png", width = 17, height = 12)

my_sens_plot(m = "FLake", l = "Stechlin", res_cali = res_cali,
                           res_sens = res_o)
ggsave("Plots/FLake_stechlin.png", width = 17, height = 12)

my_sens_plot(m = "Simstrat", l = "Erken", res_cali = res_cali,
                          res_sens = res_o)
ggsave("Plots/Simstrat_erken.png", width = 17, height = 12)


# look at some of the lakes with large interaction measure
dat_iat |> filter(iat > 0.35) |> select(lake, model, var, iat) |>
  print(n = Inf)

my_sens_plot(m = "FLake", l = "Allequash", res_cali = res_cali,
             res_sens = res_o, n_contour = 7, contour = FALSE)
ggsave("Plots/FLaket_allequash.png", width = 17, height = 12)

my_sens_plot(m = "GOTM", l = "Bosumtwi", res_cali = res_cali,
             res_sens = res_o, smet = "rmse")

my_sens_plot(m = "FLake", l = "Delavan", res_cali = res_cali,
             res_sens = res_o, smet = "nse")

my_sens_plot(m = "Simstrat", l = "Tarawera", res_cali = res_cali,
             res_sens = res_o, smet = "nse")
ggsave("Plots/Simstrat_erken.png", width = 17, height = 12)

my_sens_plot(m = "Simstrat", l = "Chao", res_cali = res_cali,
             res_sens = res_o, smet = "r")


my_sens_plot(m = "GLM", l = "Tarawera", res_cali = res_cali,
             res_sens = res_o, smet = "nse", contour = FALSE)

# look at some of the lakes with small interaction measure
dat_iat |> filter(iat > 0.35) |> select(lake, model, var, iat) |>
  print(n = Inf)

my_sens_plot(m = "GLM", l = "Zurich", res_cali = res_cali,
             res_sens = res_o, smet = "nse", contour = TRUE, n_contour = 7)

my_sens_plot(m = "Simstrat", l = "Vendyurskoe", res_cali = res_cali,
             res_sens = res_o, smet = "r", contour = TRUE, n_contour = 7)

my_sens_plot(m = "GOTM", l = "Rotorua", res_cali = res_cali,
             res_sens = res_o, smet = "rmse", contour = TRUE, n_contour = 7)

my_sens_plot(m = "FLake", l = "Tahoe", res_cali = res_cali,
             res_sens = res_o, smet = "nse", contour = TRUE, n_contour = 7)
