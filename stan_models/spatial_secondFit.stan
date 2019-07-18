data {
  int<lower = 1> N;  // number of species
  vector[N] x;       // x-axis vals
  vector[N] SD;      // observation error in y-vals
  vector[N] y;       // y_vals
}

parameters {
  real <lower = -0.1,  upper = 0.1> beta0; // y-intercept
  real <lower = -0.1,  upper = 0.1> beta1; // slope wrt x
  real <lower =  0.01, upper = 0.1> sigma; // group sd
  vector[N] yhat;
}

transformed parameters{
  real mu[N]; // community mean
  
  // calculate the community mean
  for (i in 1:N) {
    mu[i] = beta0 + beta1*x[i];
  }
}

model {
  // means of each species about population mean
  for (i in 1:N) {
    yhat[i] ~ normal(mu[i], sigma);
  }

  // add observation error
  for (i in 1:N) {
    target += normal_lpdf(y[i] | yhat[i], SD[i]);
  }
}

