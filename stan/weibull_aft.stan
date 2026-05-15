// -----------------------------------------------------------------------------
// Weibull accelerated-failure-time (AFT) model for right-censored survival.
//
// Parameterization: t_i ~ Weibull(shape, scale_i) with
//   scale_i = exp(intercept + X_i * beta)
// Under AFT, beta_k > 0 means a longer time-to-event for a 1-unit increase in
// covariate k. exp(beta) is the "time ratio" (analogue of HR's inverse).
// -----------------------------------------------------------------------------
data {
  int<lower=0> N;                          // number of subjects
  vector<lower=0>[N] times;                // observed times (event or censoring)
  array[N] int<lower=0, upper=1> status;   // 1 = event, 0 = right-censored
  int<lower=0> K;                          // number of covariates
  matrix[N, K] X;                          // covariate matrix
}
parameters {
  real intercept;                          // baseline log-scale
  vector[K] beta;                          // AFT log-time-ratio coefficients
  real<lower=0> shape;                     // Weibull shape (> 0)
}
model {
  // Weakly informative priors
  intercept ~ normal(0, 5);
  beta      ~ normal(0, 2);
  shape     ~ gamma(2, 0.5);               // E[shape] = 4; admits shape <1 or >1

  for (n in 1:N) {
    real scale = exp(intercept + X[n] * beta);
    if (status[n] == 1) {
      target += weibull_lpdf(times[n] | shape, scale);
    } else {
      target += weibull_lccdf(times[n] | shape, scale);
    }
  }
}
generated quantities {
  // exp(beta_k) is the time ratio for covariate k (>1 means survival extended)
  vector[K] time_ratio = exp(beta);
  // Posterior predictive draws for KM-overlay PP check (one t per subject per draw).
  vector[N] y_rep;
  for (n in 1:N) {
    real scale_pp = exp(intercept + X[n] * beta);
    y_rep[n] = weibull_rng(shape, scale_pp);
  }
}
