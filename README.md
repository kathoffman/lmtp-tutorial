# Introducing longitudinal modified treatment policies: a unified framework for studying complex exposures

Corresponding code guide for Hoffman et al (2023).

## Repository contents

### Main scripts: 

- `R/run_lmtp.R` - main analysis script. sets up exposure, confounders, outcome, and runs `lmtp_sdr()` for an intervention and null intervention for all time points in the study.
- `R/clean_results.R` - summarizes LMTP results. checks density ratios, and creates graphs of output.

### Supporting scripts: 

- `R/viz.R` contains data visualization functions

- 

