###############################################################################
# HIV back-calculation
# - Quarterly back-calculation with 3 CD4 bands (500+, 200-500, <200)
# - Random-walk prior on log incidence h[t] and logit diagnosis probs d[t,k]
#
# This script assumes you have already constructed the JAGS data inputs:
#   TT        : number of quarters
#   KK        : number of CD4 bands
#   q         : length-3 progression probabilities per quarter (q1,q2,q3)
#   X_HIV     : length-TT vector, non-AIDS HIV diagnoses per quarter
#   X_AIDS    : length-TT vector, AIDS diagnoses per quarter (WHO stage 4 at dx)
#   W         : TT x KK matrix, CD4-band counts among X_HIV diagnoses each quarter
#              (row sums should equal X_HIV)
###############################################################################

library(rjags)
library(coda)

# =======================
# 0) MCMC settings
# =======================
CHAINS <- 4
ITER   <- 6000
BURNIN <- 3000
THIN   <- 1

# =======================
# 1) JAGS model
# =======================
model_string <- "
model {
  # ---- initial pools and initial coverage ----
  for(k in 1:KK) { E_k0[k] ~ dgamma(1,0.001) }
  coverage_0 ~ dnorm(0.10, 100) T(0,1)          
  sumE <- sum(E_k0)
  C_0  <- sumE / (1 - coverage_0)               
  D_0  <- coverage_0 * C_0                      

  # ---- quarterly incidence random walk (log-scale) ----
  lowerG <- log(10); upperG <- log(100000)
  gamma[1] ~ dnorm(log(3200), 1/(1.0^2)) T(log(10), log(100000))
  tau_h ~ dgamma(1000, 10)
  for(i in 2:TT) {
    gamma[i] ~ dnorm(gamma[i-1], tau_h) T(lowerG, upperG)
  }
  for(i in 1:TT) { h[i] <- exp(gamma[i]) }

  # ---- diagnosis probabilities random walk (logit-scale), by CD4 band ----
  for(k in 1:KK) {
    delta[1,k] ~ dnorm(logit(0.05), 1)
    tau_d[k]   ~ dgamma(1, 0.1)
    for(i in 2:TT) {
      delta[i,k] ~ dnorm(delta[i-1,k], tau_d[k])
    }
    for(i in 1:TT) {
      d[i,k] <- ilogit(delta[i,k])
    }
  }

  # ---- cohort evolution and observation model ----
  for(k in 1:KK) { E_k[1,k] <- E_k0[k] }

  for(i in 1:TT) {
    for(k in 1:KK) { mu_HIV_k[i,k] <- E_k[i,k] * d[i,k] }
    mu_HIV[i]  <- sum(mu_HIV_k[i,])
    mu_AIDS[i] <- E_k[i,3] * (1 - d[i,3]) * q[3]

    X_HIV[i]  ~ dpois(mu_HIV[i])
    X_AIDS[i] ~ dpois(mu_AIDS[i])

    # CD4 composition among non-AIDS HIV diagnoses
    W[i,1:KK] ~ dmulti(mu_HIV_k[i,1:KK] / mu_HIV[i], X_HIV[i])

    # state transitions (undx pools)
    E_k[i+1,1] <- h[i] +
      E_k[i,1] * (1 - d[i,1]) * (1 - q[1])

    E_k[i+1,2] <- E_k[i,1] * (1 - d[i,1]) * q[1] +
      E_k[i,2] * (1 - d[i,2]) * (1 - q[2])

    E_k[i+1,3] <- E_k[i,2] * (1 - d[i,2]) * q[2] +
      E_k[i,3] * (1 - d[i,3]) * (1 - q[3])
  }

  # ---- cumulative infections/diagnoses and coverage ----
  C[1] <- C_0 + h[1]
  D[1] <- D_0 + mu_HIV[1] + mu_AIDS[1]
  coverage[1] <- D[1] / C[1]

  for(i in 2:TT) {
    C[i] <- C[i-1] + h[i]
    D[i] <- D[i-1] + mu_HIV[i] + mu_AIDS[i]
    coverage[i] <- D[i] / C[i]
  }

  # ---- annual summaries ----
  nYears <- TT %/% 4
  for(j in 1:nYears) {
    yearly_h[j]    <- sum(h[(4*j-3):(4*j)])
    yearly_diag[j] <- sum(mu_HIV[(4*j-3):(4*j)]) + sum(mu_AIDS[(4*j-3):(4*j)])
  }
}
"

# =======================
# 2) DATA: provide these objects from your own preprocessing
# =======================
# REQUIRED (define these before running):
#   TT, KK, q, X_HIV, X_AIDS, W
stopifnot(exists("TT"), exists("KK"), exists("q"), exists("X_HIV"), exists("X_AIDS"), exists("W"))

data_list <- list(
  TT    = TT,
  KK    = KK,
  q     = q,
  X_HIV = X_HIV,
  X_AIDS= X_AIDS,
  W     = W
)

# =======================
# 3) Compile and sample
# =======================
jmod <- jags.model(
  textConnection(model_string),
  data     = data_list,
  n.chains = CHAINS,
  n.adapt  = 2000,
  quiet    = TRUE
)

update(jmod, BURNIN)

samps <- coda.samples(
  jmod,
  variable.names = c(
    "mu_HIV","mu_AIDS","mu_HIV_k","d",
    "C","D","coverage",
    "yearly_h","yearly_diag",
    "C_0","D_0",
    "tau_h","tau_d",
    "gamma","delta","h"
  ),
  n.iter = ITER,
  thin   = THIN
)

