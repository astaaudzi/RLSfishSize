data {
  int<lower = 1> N;  // number of data
  vector[N] dd;      // decimal date (since 2018)
  vector[N] msst; // mean sea-surface temperature
}

parameters {
  real <lower = 10.0,  upper = 30.0> beta0; // mean temperature
  real <lower = -2.0,  upper =  2.0> beta1; // annual rate of change
  real <lower = 0.01, upper =  10.0> beta2; // annual amplitude
  real <lower = -0.40, upper =  0.6> phi;   // peak time
  real <lower =  0.01, upper = 10.0> sigma; // observation error
}

model {
  msst ~ normal(beta0 + beta1 * dd + beta2 * cos(2*3.14159*(dd-phi)), sigma);
}
