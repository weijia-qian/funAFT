## code to prepare `pupil` dataset goes here
load("data-raw/df_post1_r.rda")

pupil <- df_post1_r[, !names(df_post1_r) %in% c("trial_id", "pct_chg_120")]

usethis::use_data(pupil, overwrite = TRUE)
