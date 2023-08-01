library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)


res <- read.csv("res_sens.csv")




res |> filter(lake == "Annie" & var == "rmse") |> ggplot() +
  geom_col(aes(y = delta, x = names)) +
  geom_errorbar(aes(x = names, ymin = delta - delta_conf/2,
                    ymax = delta + delta_conf), col = 2) +
  theme_pubr() + grids() + facet_wrap(~model, scales = "free_x") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  xlab("parameter") + ggtitle("Annie")

res_mip <- res |> group_by(lake, model) |>
  reframe(par_d = names[delta == max(delta)],
          par_S1 = names[S1 == max(S1)])

res_mip |> pivot_longer(cols = 3:4) |> ggplot() +
  geom_histogram(aes(x = value, fill = name), stat = "count", position = "dodge") +
  facet_wrap(~model, scales = "free_x") + theme_pubr() + grids() +
  xlab("parameter")

# check how many times delta and S1 differ
sum(res_mip$par_d != res_mip$par_S1)
res_mip[res_mip$par_d != res_mip$par_S1, ]
