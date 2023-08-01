library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)


res <- read.csv("results/res_sens.csv")




res |> filter(lake == "Annie" & var == "rmse") |> ggplot() +
  geom_col(aes(y = delta, x = names)) +
  geom_errorbar(aes(x = names, ymin = delta - delta_conf/2,
                    ymax = delta + delta_conf), col = 2) +
  theme_pubr() + grids() + facet_wrap(~model, scales = "free_x") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  xlab("parameter") + ggtitle("Annie")

res |> group_by(lake, model) |> summarise(par_s = names[delta == max(delta)]) |>
  ggplot() + geom_histogram(aes(x = par_s), stat = "count") +
  facet_wrap(~model, scales = "free_x") + theme_pubr() + grids() +
  xlab("parameter")
