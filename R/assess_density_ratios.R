# density ratios

library(lmtp) 
library(tidyverse)
library(survival)
library(gt)
library(wesanderson)
library(viridis)

set.seed(7)

dat <- read_rds(here::here("data/toy_data.rds"))
out_shift <- read_rds(here::here("results/out_shift_14.rds"))
out_null <- read_rds(here::here("results/out_null_14.rds"))

# find high density ratios ------------------------------------------------------

get_drs <- function(result_file, data_file){
    dr <- as.data.frame(result_file[["density_ratios"]]) 
    colnames(dr) <- paste0("density_ratio_",1:ncol(dr))
    dr_merg <- bind_cols(data_file, dr) |>
        mutate(id = factor(row_number()))
    return(dr_merg)
}

extract_high_drs <- function(result_file, data_file, cutoff){
    get_drs(result_file, data_file) %>%
        select(id, starts_with("density_ratio")) %>%
        pivot_longer( starts_with("density_ratio")) %>%
        group_by(id) %>%
        mutate(max_density_ratio = max(value)) %>%
        mutate(time = factor(parse_number(name))) %>%
        filter(max_density_ratio > cutoff) 
}

plot_high_drs <- function(result_file, data_file, cutoff){
    extract_high_drs(result_file, data_file, cutoff) %>%
        ggplot(aes(time, value)) +
        geom_jitter() +
        geom_violin() +
        geom_boxplot(width=.2, fill=NA, outlier.colour = NA) +
        theme_classic() +
        labs(x = "time", y="density_ratio")
}

# plot density ratios across time ------------------------------------------------------

# note that a density ratio of exactly 0 indicates that observation was censored,
# so we could filter these out
# note also that a density ratio close to 0 indicates a similar prob of receiving
# the intervened exposure, conditional on covariates, as receiving the naturally
# observed exposure, conditional on covariates. not a positivity violation (unlike
# a near zero propensity score

get_drs(out_shift, dat)  %>%
    select(id, starts_with("density_ratio")) %>%
    pivot_longer( starts_with("density_ratio")) %>% # filter(row_number() == 3) %>% select(value)
    filter(value > 0) %>%
    mutate(time = parse_number(name)) %>%
    ggplot(aes(value)) +
    geom_histogram() +
    facet_wrap(~time, scales="free") +
    labs(title = "Distribution of raw density ratios") +
    theme_bw()
