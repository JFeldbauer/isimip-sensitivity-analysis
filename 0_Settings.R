# settings for the scripts to analyse the ISIMIP LER calibration runs

# plotting theme
thm <- theme_pubr(base_size = 16) + grids()

# which performance measures to include
p_metrics <- c("bias",
               "mae",
               #"nmae",
               "nse",
               "r",
               "rmse")
