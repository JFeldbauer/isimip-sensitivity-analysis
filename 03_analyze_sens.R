setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# clean up
rm(list=ls())
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

# set sensitivity values that are smaller than the of the dummy variable
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

##--------- single most sensitive parameter  ---------------------------

# calculate single most sensitive parameter for each lake and model
res_mip <- res |> group_by(lake, model, var) |>
  reframe(par_d = names[delta == max(delta)],
          par_S1 = names[S1 == max(S1)])

# plot distribution of single most sensitive parameter over all lakes and models
# for all measures
res_mip |> pivot_longer(cols = 4:5) |>
  mutate(name = case_match(name,
                           "par_d" ~ "\u03B4",
                           "par_S1" ~ "S1")) |> ggplot() +
  geom_histogram(aes(x = value, fill = name), stat = "count", position = "dodge") +
  facet_grid(var~model, scales = "free_x") + thm +
  grids() + xlab("parameter") +
  theme(axis.text.x=element_text(angle = -60, hjust = 0)) +
  scale_fill_manual("Sensitivity \n measure", values = c("#45B2DD", "#72035F"))

ggsave("Plots/sens_single.png", width = 11, height = 9)

# only use delta sensitivity metric but look at cluster
res_mip |> pivot_longer(cols = par_d) |>
  mutate(name = case_match(name,
                           "par_d" ~ "\u03B4",
                           "par_S1" ~ "S1")) |>
  left_join(meta, by = c("lake" = "Lake.Short.Name")) |> ggplot() +
  geom_histogram(aes(x = value, fill = kmcluster),
                 stat = "count", position = "dodge", col = 1) +
  facet_grid(var~model, scales = "free_x") + thm +
  grids() + xlab("parameter") +
  theme(axis.text.x=element_text(angle = -60, hjust = 0)) +
  scale_fill_viridis_d("Cluster")


# check if there is a correlation between sensitive parameter and cluster
res_mip |> left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  select(model, var, par_d, par_S1, kmcluster) |>
  pivot_longer(3:4) |>
  mutate(name = case_match(name,
                           "par_d" ~ "\u03B4",
                           "par_S1" ~ "S1")) |> ggplot() +
  geom_histogram(aes(x = value, fill = name), stat = "count", position = "dodge") +
  facet_grid(paste("Cluster ",kmcluster)~model, scales = "free") +
  theme_pubr(base_size = 17) +
  theme(axis.text.x=element_text(angle = -55, hjust = 0)) +
  grids() + xlab("parameter") +
  scale_fill_manual("Sensitivity \n measure", values = c("#45B2DD", "#72035F"))

ggsave("Plots/sens_single_clust.png", width = 11, height = 9)

# check how many times most sensitive parameter from delta and S1 differ
sum(res_mip$par_d != res_mip$par_S1)
res_mip[res_mip$par_d != res_mip$par_S1, ]

## relate most sensitive parameter to meta info on lake
res_mip <- left_join(res_mip, meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(Reservoir.or.lake. = as.factor(Reservoir.or.lake.),
         par_d = as.factor(par_d),
         par_S1 = as.factor(par_S1))

res_mip |> ggplot() + geom_point(aes(y = par_d, x = mean.depth.m,
                                     col = var), size = 3) +
  facet_wrap(~model, scales = "free_y") + scale_x_log10() + thm

##--------- group of most sensitive parameters -------------------------------

# calculate group of most sensitive parameters for each lake and model
# function to return all values in a vector that contribute to frac (default 75)
# percent of the sum of all values
gmiv <- function(x, frac = .75) {
  tmp <- sort(x, decreasing = TRUE)

  
  id <- which(cumsum(tmp)/sum(tmp) <= frac)
  if(length(id) > 0) {
    im <- (max(id)) + 1
  } else {
    im <- 1
  }

  id <- 1:im
  return(tmp[id])
}

# for each model, lake, and measure return the most sensitive parameters
res_gip <- res |> group_by(lake, model, var) |> filter(delta != 0 & S1 != 0) |>
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
delta_gip <- res_gip |> pivot_longer(seq(4, 14, by = 2)) |>
  filter(grepl("par_d_.*", name)) |>
  mutate(frac = gsub(".*_", "", name)) |> 
  group_by(model, var, frac, lake) |> reframe(par = unlist(strsplit(value, ", ")),
                                                   meas = "\u03B4")
S1_gip <- res_gip |>pivot_longer(seq(4, 14, by = 2)) |>
  filter(grepl("par_S1_.*", name)) |>
  mutate(frac = gsub(".*_", "", name)) |> 
  group_by(model, var, frac, lake) |> reframe(par = unlist(strsplit(value, ", ")),
                                                meas = "S1")
rbind(delta_gip, S1_gip) |> filter(frac == "1") |>
  ggplot() +
  geom_histogram(aes(x = par, fill = meas), stat = "count", position = "dodge") +
  facet_grid(var~model, scales = "free_x") + thm +
  grids() + xlab("parameter") +
  theme(axis.text.x=element_text(angle = -55, hjust = 0)) +
  scale_fill_manual("Sensitivity \n measure", values = c("#45B2DD", "#72035F"))

ggsave("Plots/count_sens.png", width = 14, height = 11)

# also plot for the different cluster
rbind(delta_gip, S1_gip) |> filter(frac == "1") |>
  left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  ggplot() + geom_histogram(aes(x = par,  fill = meas),
                            stat = "count", position = "dodge") + 
  facet_grid(kmcluster~model, scales = "free_x") + thm + grids() +
  scale_fill_manual("Sensitivity \n measure", values = c("#45B2DD", "#72035F")) +
  theme(axis.text.x=element_text(angle = -55, hjust = 0)) +
  xlab("Parameter")

ggsave("Plots/count_sens_clust.png", width = 14, height = 11)

# same but sum up the different performacne metrics
rbind(delta_gip, S1_gip) |> filter(frac == "1") |>
  left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  group_by(model, kmcluster, var, par) |> reframe(cnt = length(var)) |>
  group_by(model, kmcluster, var) |> mutate(cnt = cnt/sum(cnt)) |>
  ungroup() |>
  #complete(model, kmcluster, var, par) |>
  ggplot() + geom_tile(aes(x = par,  y = kmcluster, fill = cnt)) + 
  facet_grid(var~model, scales = "free_x") + thm + grids() +
  scale_fill_viridis_c("Frequency", option = "C") +
  theme(axis.text.x=element_text(angle = -55, hjust = 0)) +
  xlab("Parameter") + ylab("Cluster")

ggsave("Plots/freq_sens_clust.png", width = 14, height = 11)


## plot distribution of number of sensitive parameters
rbind(delta_gip, S1_gip) |> group_by(model, var, frac, lake, meas) |>
  reframe(n = n()) |> mutate(frac = factor(frac, levels = c("1", "75", "05"),
                                           labels = c("100%", "75%", "50%"))) |>
  ggplot() + geom_histogram(aes(x = n,  fill = frac),
                            stat = "count", position = "dodge") + 
  facet_grid(var~model) + thm + grids() +
  scale_fill_viridis_d("Fraction of total sum", option = "D") +
  xlab("Number of parameters contributing")

ggsave("Plots/count_imp_par.png", width = 14, height = 11)

## alternative plot
rbind(delta_gip, S1_gip) |> group_by(model, var, frac, lake, meas) |>
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

# lakes where no parameter is sensitive
res |> group_by(lake, var, model) |>
  mutate(no_sens = all(delta[names != "dummy"] <= delta[names == "dummy"]) ||
           all(S1[names != "dummy"] <= S1[names == "dummy"])) |>
  filter(no_sens) |> select(lake) |> unique()  # only happens for nmae


# only for nmae there are cases where no parameter is sensitive


##------------------- look at sensitivity of wind speed scaling ---------------

res |> left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(delta = ifelse(delta > 1, 0, delta)) |>
  filter(names == "wind_speed") |> ggplot() +
  geom_point(aes(x = lake.area.sqkm, y = mean.depth.m, col = delta)) +
  facet_grid(model~var) +
  scale_x_log10() + scale_y_log10() + scale_color_viridis_c(option = "D")

res |> left_join(meta, by = c("lake" = "Lake.Short.Name")) |>
  mutate(delta = ifelse(delta > 1, 0, delta)) |>
  filter(names == "wind_speed") |> ggplot() +
  geom_boxplot(aes(x = kmcluster, y = delta, fill = kmcluster)) +
  facet_grid(model~var) + thm +
  scale_fill_viridis_d() + xlab("delta sensitivity wind speed scaling")


##---------- plots for single models -----------------------------

#visual comparison of heatmaps against the calculated sensitivity metrics

# function to plot heatmaps for model performance along the different parameter
my_sens_plot <- function(m = "GLM", l = "Zurich", res_cali, res_sens,
                         smet = "rmse") {
  
  sens <- res_sens |> filter(lake == l & model == m & var == smet)
  pars <- unique(sens$names)
  pars_s <- pars
  pars <- pars[pars!="dummy"]
  log <- rep(FALSE, length(pars))
  log[pars == "turb_param.k_min"] <- TRUE
  
  spec <- viridis::viridis(15)
  dat <- filter(res_cali, model == m & lake == l)
  dat_best <- filter(dat, rmse < quantile(rmse, 0.1, na.rm = TRUE))
  thm <- theme_pubr(base_size = 13) + grids()
  
  pl <- lapply(combn(pars, 2, simplify = FALSE), function(p) {
    plt <- ggplot(dat_best) +
      geom_point(aes_string(x = p[1], y = p[2], color = smet), shape = 15,
                 size = 2, alpha = 0) +
      geom_point(data = dat,
                 aes_string(x = p[1], y = p[2], color = smet), shape = 15,
                 size = 1.5, alpha = 0.75) +
      scale_colour_gradientn(colours = rev(spec)) + thm +
      theme(legend.position = "none",
            plot.margin = margin(t = 20,    # Top margin
                                 r = 30,    # Right margin
                                 b = 10,    # Bottom margin
                                 l = 10))
      
    
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
  
  do.call(grid.arrange, c(pl2, ncol = length(pars)-1, as.table = FALSE) )
  
  
  
}

p_gtm_kivu <- my_sens_plot(m = "GOTM", l = "Kivu", res_cali = res_cali,
                           res_sens = res_o)
p_glm_biel <- my_sens_plot(m = "GLM", l = "Biel", res_cali = res_cali,
                           res_sens = res_o)
p_fl_stech <- my_sens_plot(m = "FLake", l = "Stechlin", res_cali = res_cali,
                           res_sens = res_o)
p_sim_erk <- my_sens_plot(m = "Simstrat", l = "Erken", res_cali = res_cali,
                          res_sens = res_o)



ggsave("Plots/GOTM_kivu.png", plot = p_gtm_kivu, width = 17, height = 12)
ggsave("Plots/GLM_biel.png", plot = p_glm_biel, width = 17, height = 12)
ggsave("Plots/FLake_stechlin.png", plot = p_fl_stech, width = 17, height = 12)
ggsave("Plots/Simstrat_erken.png", plot = p_sim_erk, width = 17, height = 12)
