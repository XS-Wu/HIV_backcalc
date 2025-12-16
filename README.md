# HIV Back-calculation (Core JAGS Model Code)

This repository contains core model code for a HIV back-calculation model with 3 CD4 strata at diagnosis.

✅ Included:
- JAGS model specification (`model_string`)
- JAGS compilation + MCMC sampling code (rjags + coda)
- Expected input objects and their required shapes

## Model overview

- Time index: quarters `t = 1..TT`
- CD4 strata: `k = 1..3` representing `500+`, `200–500`, `<200`
- Latent undiagnosed pools `E_k[t,k]` evolve by:
  - New infections entering the highest CD4 stratum
  - CD4 progression with fixed per-quarter probabilities `q`
  - Diagnosis with time-varying probabilities `d[t,k]`

### Priors
- Incidence `h[t]` follows a **random walk** on log scale.
- Diagnosis probabilities `d[t,k]` follow **random walks** on logit scale.

### Likelihood
- Non-AIDS HIV diagnoses `X_HIV[t] ~ Poisson(mu_HIV[t])`
- AIDS diagnoses `X_AIDS[t] ~ Poisson(mu_AIDS[t])`
- CD4-band composition among non-AIDS HIV diagnoses:
  `W[t,] ~ Multinomial(proportions[t,], X_HIV[t])`

## Requirements

### R packages
- rjags
- coda

Also requires a working JAGS installation.

## Inputs you must prepare

Before sourcing/running the core script, you must define:

- `TT` (integer): number of quarters
- `KK` (integer): number of CD4 strata
- `q` (numeric length 3): quarterly progression probs `(q1,q2,q3)`
- `X_HIV` (numeric/integer length TT): non-AIDS HIV diagnoses per quarter
- `X_AIDS` (numeric/integer length TT): AIDS diagnoses per quarter
- `W` (matrix TT x KK): CD4-band counts among `X_HIV` each quarter  

## How to run

1. Prepare the six objects above (`TT`, `KK`, `q`, `X_HIV`, `X_AIDS`, `W`).
2. Source the script:
   ```r
   source("backcalc_core.R")
   ```
3. The MCMC samples are stored in `samps` (a `coda` mcmc.list).

Optionally save:
```r
saveRDS(samps, file = "backcalc_core.rds")
```
