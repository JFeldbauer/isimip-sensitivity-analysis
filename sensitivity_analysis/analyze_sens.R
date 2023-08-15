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

# load results of sensitivity analysis
res <- read.csv("res_sens.csv")
# load results of latin hypercube calibration
res_cali <- read.csv("../data/results_lhc.csv")
# load lake meta data
meta <- read.csv("../data/Lake_meta.csv")

##----------------- first plots ---------------------------------

# plot delta sensitivity metrics for one lake
res |> filter(lake == "Annie" & var == "rmse") |> ggplot() +
  geom_col(aes(y = delta, x = names)) +
  geom_errorbar(aes(x = names, ymin = delta - delta_conf/2,
                    ymax = delta + delta_conf), col = 2) +
  theme_pubr() + grids() + facet_wrap(~model, scales = "free_x") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  xlab("parameter") + ggtitle("Annie")


##--------- single most sensitive parameter  ---------------------------

# calculate single most sensitive parameter for each lake and model
res_mip <- res |> group_by(lake, model, var) |>
  reframe(par_d = names[delta == max(delta)],
          par_S1 = names[S1 == max(S1)])

# plot distribution of single most senstitive parameter over all lakes and models
# for all measures
res_mip |> pivot_longer(cols = 4:5) |>
  mutate(name = case_match(name,
                           "par_d" ~ "delta",
                           "par_S1" ~ "S1")) |> ggplot() +
  geom_histogram(aes(x = value, fill = name), stat = "count", position = "dodge") +
  facet_wrap(~model, scales = "free_x") + theme_pubr(base_size = 17) +
  grids() + xlab("parameter") +
  scale_fill_manual("Sensitivity \n measure", values = c("#45B2DD", "#72035F"))

ggsave("count_sens.png", width = 14, height = 9)

# check how many times most sensitive parameter from delta and S1 differ
sum(res_mip$par_d != res_mip$par_S1)
res_mip[res_mip$par_d != res_mip$par_S1, ]

## relate most sensitive parameter to meta info on lake
res_mip <- left_join(res_mip, meta, by = c("lake" = "Lake.Short.Name")) |>
  select(-Country, -Lake.Name.Folder, -Duration, -Lake.Name) |>
  mutate(Reservoir.or.lake. = as.factor(Reservoir.or.lake.),
         par_d = as.factor(par_d),
         par_S1 = as.factor(par_S1))

res_mip |> ggplot() + geom_point(aes(y = par_d, x = mean.depth.m,
                                     col = var), size = 3) +
  facet_wrap(~model, scales = "free_y") + scale_x_log10()

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
res_gip <- res |> group_by(lake, model, var) |>
  reframe(par_d = paste0(names[delta %in% gmiv(delta)], collapse = ", "),
          n_par_d = length(gmiv(delta)),
          par_S1 = paste0(names[S1 %in% gmiv(S1)], collapse = ", "),
          n_par_S1 = length(gmiv(S1)))

# plot distribution of most senstitive parameters over all lakes and models
# for both measures
delta_gip <- res_gip |> group_by(model, var) |> reframe(par = unlist(strsplit(par_d, ", ")),
                                                   meas = "delta")
S1_gip <- res_gip |> group_by(model, var) |> reframe(par = unlist(strsplit(par_S1, ", ")),
                                                meas = "S1")
rbind(delta_gip, S1_gip) |>
  ggplot() +
  geom_histogram(aes(x = par, fill = meas), stat = "count", position = "dodge") +
  facet_grid(var~model, scales = "free_x") + theme_pubr(base_size = 17) +
  grids() + xlab("parameter") +
  theme(axis.text.x=element_text(angle = -55, hjust = 0)) +
  scale_fill_manual("Sensitivity \n measure", values = c("#45B2DD", "#72035F"))

ggsave("count_sens.png", width = 14, height = 9)



##---------- plots for single models -----------------------------

#visual comparison of heatmaps against the calculated sensitivity metrics

# function to plot heatmaps for model performance along the different parameter
my_sens_plot <- function(m = "GLM", l = "Zurich", res_cali, res_sens,
                         smet = "rmse") {
  
  sens <- res_sens |> filter(lake == l & model == m & var == smet)
  pars <- unique(sens$names)
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
  
  ps1 <- sens |> ggplot() +
    geom_col(aes(y = delta, x = names), fill = "#45B2DD") +
    geom_errorbar(aes(x = names, ymin = delta - delta_conf/2,
                      ymax = delta + delta_conf),
                  col = "#0D3D20", lwd = 1.25, width=.5) +
    theme_pubr() + grids() +
    theme(axis.text.x=element_text(angle = -65,
                                   hjust = 0, size = 9),
          plot.margin = margin(t = -40,    # Top margin
                               r = 10,    # Right margin
                               b = -70,    # Bottom margin
                               l = 10)) + # Left margin
    xlab("") + ylab("")
    
  ps2 <- sens |> ggplot() +
      geom_col(aes(y = S1, x = names), fill = "#72035F") +
      geom_errorbar(aes(x = names, ymin = S1 - S1_conf/2,
                        ymax = S1 + S1_conf),
                    col = "#03DE67", lwd = 1.25, width=.5) +
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

  pl2[[21]] <- as_ggplot(text_grob("Delta"))
  pl2[[16]] <- as_ggplot(text_grob("S1"))
  pl2[[22]] <- ps1
  pl2[[17]] <- ps2
  pl2[[6]] <- t
  
  do.call(grid.arrange, c(pl2, ncol = length(pars)-1, as.table = FALSE) )
  
  
  
}

p_gtm_kivu <- my_sens_plot(m = "GOTM", l = "Kivu", res_cali = res_cali, res_sens = res)
p_glm_biel <- my_sens_plot(m = "GLM", l = "Biel", res_cali = res_cali, res_sens = res)
p_fl_stech <- my_sens_plot(m = "FLake", l = "Stechlin", res_cali = res_cali, res_sens = res)
p_sim_erk <- my_sens_plot(m = "Simstrat", l = "Erken", res_cali = res_cali, res_sens = res)



ggsave("GOTM_kivu.png", plot = p_gtm_kivu, width = 17, height = 12)
ggsave("GLM_biel.png", plot = p_glm_biel, width = 17, height = 12)
ggsave("FLake_stechlin.png", plot = p_fl_stech, width = 17, height = 12)
ggsave("Simstrat_erken.png", plot = p_sim_erk, width = 17, height = 12)
