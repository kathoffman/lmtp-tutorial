## ------------------------------------------------------------------------------------------
##
## MAIN SCRIPT: RUN LMTP
##
## Purpose of script: Show an example of code to run an analysis for a 
##     Modified Treatment Policy (MTP) for a time-varying, multi-level exposure,
##          time-dependent confounders, informative right censoring for a time-to-event outcome.
##
## Author: Katherine Hoffman
##
## Date Created: 2023-04-11
##
## Author: Katherine Hoffman, 2023
## Email: kathoffman.stats@gmail.com
##
## ------------------------------------------------------------------------------------------
##
## Notes: The data set of n=52000 patients for time=(1,2,...,14) contains realistic values for patients but
## is a toy data set and cannot be used to obtain the result in the paper due to
## data-sharing restrictions. The interventions are
##.     1: A one-day delay in intubation, no loss to follow up
##      2: Intubation as observed, no loss to follow up
##
## ------------------------------------------------------------------------------------------

# load libraries, set seed, read in functions ---------------------------------------------------

library(lmtp)
library(tidyverse)
library(future)
library(earth)
library(ranger)
library(xgboost)

set.seed(7)

dat <- read_rds("data/toy_data.rds") 

# set up parameters to input to lmtp_sdr function ---------------------------------------------------

# for calling column names
outcome_day <- 14
padded_days <- str_pad(0:(outcome_day - 1), 2, pad = "0")
padded_days_out <- str_pad(1:outcome_day, 2, pad = "0")

# set parameters of outcome, trt, censoring and baseline adjustment vars
a <-  paste0("I_", padded_days) # treatment
bs <- dat |> dplyr::select(age:active_cancer) |> names()
y <- paste0("Y_",padded_days_out) # outcome (AKI)
censoring <- paste0("C_",padded_days) # observed at next time

used_letters <- dat |> # get all letters for time varying
  dplyr::select(starts_with("L_")) |>
  names()

# turn time varying col names into a list (see ?lmtp_sdr tv argument input)
tv <- map(0:(outcome_day - 1), function(x) { # list for time varying covars
  used_letters[str_detect(used_letters, str_pad(x, 2, pad = "0"))]
})

# write custom shift function ---------------------------------------------------

# A modified treatment policy (MTP) to delay everyone's intubation (I_* == 2) by 1 day.
# Instead, set them to non-invasive supp O2 (I_* == 1)
mtp <- function(data, trt) {
  # extract time point
  tau <- readr::parse_number(trt)
  # get the col name of previous trt
  trt_prev <- paste0("I_", stringr::str_pad(tau - 1, 2, "left", "0"))
  
  if(trt == "I_00") {
    # if first time point and intubated, set to 1
    data[data[[trt]] == "2", trt] <- factor("1", levels=0:2)
  } else {
    # if intubated at time T but not T-1, set to 1
    data[which(data[[trt]] == "2" & data[[trt_prev]] != "2"), trt] <- factor("1", levels=0:2)
  }
  return(factor(data[[trt]], levels = 0:2)) # return the refactored treatment level
}

# set up analysis parameters ----------------------------------------------

# useful constants
trim <- 0.995     # density ratio trimming - akin to truncating a propensity score
folds <- 5        # "outer" folds for cross-fitting
SL_folds <- 5     # "inner" folds for super learning
k <- 2            # how much history is used at each t (markov assumption)

lrnrs <-  c("SL.glm", "SL.earth", "SL.xgboost", "SL.ranger") # superlearner candidate libraries

#################################################################################
# Run LMTP functions for 1. shift (delay in intubation + no loss to follow up)
# and 2. no shift (intubation as observed, no loss to follow up)
# for all time points (outcome = day 1, 2, ... 14)
# if we just want survival differences for day 14, only need to run for day 14
# but we wanted to estimate the entire survival curve
#################################################################################

for (this_time in 14:1) { # could be switched to foreach, future_map, etc.

progressr::with_progress(
out_shift <-
  lmtp_sdr(
    dat,
    trt = a[seq_len(this_time)],
    outcome = y[seq_len(this_time)],
    baseline = bs,
    time_vary = tv[seq_len(this_time)],
    cens = censoring[seq_len(this_time)],
    shift = mtp, # call my custom mtp function to delay intubation
    mtp = T, # tell LMTP package this is an MTP
    outcome_type = ifelse(this_time == 1, "binomial", "survival"),
    learners_outcome = lrnrs,
    learners_trt = lrnrs,
    folds = folds,
    .SL_folds = SL_folds,
    .trim = trim,
    k = k
  )
)

# save the out file for the MTP shift
saveRDS(out_shift, paste0("results/out_shift_", this_time, ".rds"))

progressr::with_progress(
  out_null <-
    lmtp_sdr(
      dat,
      trt = a[seq_len(this_time)],
      outcome = y[seq_len(this_time)],
      baseline = bs,
      time_vary = tv[seq_len(this_time)],
      cens = censoring[seq_len(this_time)],
      outcome_type = ifelse(this_time == 1, "binomial", "survival"),
      learners_outcome = lrnrs,
      learners_trt = lrnrs,
      folds = folds,
      .SL_folds = SL_folds,
      .trim = trim,
      k = k
      )
)

# save the out file for the null shift
saveRDS(out_null, paste0("results/out_null_", this_time, ".rds"))
        
}

