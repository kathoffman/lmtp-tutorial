## ------------------------------------------------------------------------------------------
##
## MAIN SCRIPT: CLEAN LMTP RESULTS
##
## Purpose of script: Clean output (into tables and plots) from LMTP objects
## 
## Author: Katherine Hoffman
##
## Date Created: 2023-04-11
##
## Author: Katherine Hoffman, 2023
## Email: kathoffman.stats@gmail.com
##
## ------------------------------------------------------------------------------------------

# load libraries, setup -------------------------------------------------------------

library(tidyverse)
library(lmtp)

source(here::here("R/utils.R")) # source helper functions
source(here::here("R/vis.R"))

# read in results (lmtp objects) --------------------------------------------------

# Read in LMTP results for MTP and null shift as two lists
est_intubation <- list()
est_control <- list()
for (i in 1:14){
  tmp <- read_rds(here::here(paste0("results/out_shift_", i, ".rds")))
  est_intubation[[i]] <- tmp
  tmp <- read_rds(here::here(paste0("results/out_null_", i, ".rds")))
  est_control[[i]] <- tmp
}

# switch the first time point (t=1) from incidence to survival estimates
est_intubation[[1]] <- switch_binary(est_intubation[[1]])
est_control[[1]] <- switch_binary(est_control[[1]])

# results as tables (tibbles) --------------------------------------------------------

# clean up results using helper functions from utils.R
# compute marginal or simulaneous confidence bands (in paper, we report simult)
marginal_results <- summarize_results(est_trt = est_intubation, est_ctl = est_control)
simult_results <- summarize_results(est_trt = est_intubation, est_ctl = est_control, ci_type = "simult")

# plotting results --------------------------------------------------------

# create graphics of results using functions in vis.R

# A) incidence estimates across time
p_surv_sdr <- simult_results$surv_est |>
  bind_rows(.id = "trt_type") |>
  mutate(
    trt_type = case_when(
      trt_type == "1" ~ "Delayed intubation (MTP)",
      trt_type == "2" ~ "No intervention"
    )
  ) |>
  plot_surv()

ggsave(p_surv_sdr, width = 12, height = 8,
       file = here::here("graphs", "sdr_surv_est_mort.pdf"))

# B) difference in incidence estimates between shifts across time
p_survdiff_sdr <- simult_results$diff_est |>
  mutate(
    p_adj = p.adjust(pval, "bonferroni"),
  ) |>
  dplyr::select(-std_err, -test_stat, -pval) |>
  plot_survdiff()

ggsave(p_survdiff_sdr, width = 12, height = 8,
       file = here::here("graphs", "sdr_survdiff_est_mort.pdf"))

# paneled plot (A, B)
p_surv_sdr_paneled <- p_surv_sdr + p_survdiff_sdr +
  plot_layout(ncol=2) + plot_annotation(tag_levels = "A")

ggsave(p_surv_sdr_paneled, width = 20, height = 9,
       file = here::here("graphs", "sdr_surv_paneled_mort.pdf"))


