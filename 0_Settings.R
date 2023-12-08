# settings for the scripts to analyse the ISIMIP LER calibration runs

# plotting theme
thm <- theme_pubr(base_size = 16) + grids()

# which performance measures to include
p_metrics <- c("bias",
               #"mae",
               #"nmae",
               "nse",
               "r",
               "rmse")

# number of clusters to create
nclust <- 5

# names of the parameter that were calibrated
par_names <- c("wind_speed", "swr", "Kw", "c_relax_C", "fetch_lk",
               "depth_bs_lk", "mixing.coef_mix_hyp", "mixing.coef_mix_conv",
               "mixing.coef_mix_turb", "turb_param.k_min", "bottom.h0b",
               "turb_param.const_num", "a_seiche", "hgeo", "cd")


# included models
models_used <- c("FLake", "GLM", "GOTM", "Simstrat")
