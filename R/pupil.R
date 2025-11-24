#' Pupil Light Response Dataset
#'
#' A dataset containing pupil light–response curves and cannabis user group information
#' for 123 subjects. Each row represents one subject.
#'
#' @details
#' The dataset includes:
#'
#' - `ptid`:
#'   Patient ID.
#'
#' - `eye`:
#'   The eye being tested (`Right`).
#'
#' - `time_since_use`:
#'   Survival outcome. Time since most recent cannabis use (in minutes).
#'
#' - `is_user`:
#'   Event indicator. `1 = observed user`, `0 = censored`.
#'
#' - `pct_chg_1` to `pct_chg_119`:
#'   Percent change in right-eye pupil diameter at 119 uniformly spaced time points
#'   during a light-stimulus test. Each row represents the functional trajectory
#'   \( X_i(s) \).
#'
#' - `age_in_years`:
#'   Subject age (years).
#'
#' - `bmi`:
#'   Body mass index.
#'
#' @format A data frame with `r nrow(pupil)` rows and `r ncol(pupil)` variables.
#'
#' @usage data(pupil)
#'
#' @source Internal dataset used in funAFT examples.
#'
#' @examples
#' data(pupil)
#' str(pupil)
#'
"pupil"
