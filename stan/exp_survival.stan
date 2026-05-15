// Exponential survival with treatment indicator. Used at adaptive interim/RAR
// updates to compute posterior on HR = lambda_treat / lambda_ctrl.
data {
  int<lower=0> N;
  vector<lower=0>[N] time;                  // observed time (event or censoring)
  array[N] int<lower=0, upper=1> evt;       // 1 = event, 0 = right-censored
  array[N] int<lower=0, upper=1> treat;     // 1 = treatment, 0 = control
}
parameters {
  real<lower=0> lambda_c;                   // control hazard per month
  real log_hr;                              // log hazard ratio (treatment vs control)
}
transformed parameters {
  real<lower=0> hr = exp(log_hr);
}
model {
  // Weakly informative priors. Prior on lambda_c is centered at the
  // baseline-event-rate truth (0.025/month = 0.30/year) so that early
  // interim posteriors with few events are not pulled toward an
  // unrealistically high hazard.
  lambda_c ~ gamma(2, 80);                  // E[lambda_c] = 2/80 = 0.025/month
  log_hr   ~ normal(0, 1);                  // 95% prior interval ~ HR in (0.14, 7.0)

  for (i in 1:N) {
    real lam = treat[i] == 1 ? lambda_c * hr : lambda_c;
    if (evt[i] == 1) {
      target += exponential_lpdf(time[i] | lam);
    } else {
      target += exponential_lccdf(time[i] | lam);
    }
  }
}
