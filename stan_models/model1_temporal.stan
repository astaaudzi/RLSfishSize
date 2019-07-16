// Shane A. Richards 20/12/2018
// fits locations independently
// includes random effects associated with year, and survey
// combines identical observations (i.e. uses counts, not single observations)

data {
  int  <lower = 1>            N;           // number of surveys
  int  <lower = 1>            J;           // number of locations
  int  <lower = 1>            K;           // number of size classes
  int  <lower = 1>            L;           // number of year classes
  real                        cutoff[K-1]; // size class cut-offs
  matrix[N,K]                 y;           // observed fish size classes (counts)
  real                        x[N];        // centred year (differs per location)
  int  <lower = 0, upper = 1> mpa[N];      // observed predictor mpa
  int  <lower = 1, upper = J> i_loc[N];    // location (index only)
  int  <lower = 1, upper = L> i_yr[N];     // observed year (index only)
}

parameters {
  real <lower =  1.0, upper = 4.0> beta_0[J];   // log(mean size) for each location
  real <lower =  -0.10, upper = 0.10> beta_yr[J];  // annual change on log(size) for each location
  real <lower =  -0.50, upper = 0.50> beta_mpa;    // change in log(mean size) when in mpa
  real <lower =  0.100, upper = 0.40> sigma_size;  // variation in (log)size
  real <lower =  0.001, upper = 0.50> sigma_yr;    // random variation among years
  real <lower =  0.001, upper = 0.50> sigma_srv;   // random variation among surveys
  real yr_RE[L];    // estimated year-specific variation (random effect)
  real srv_RE[N];   // estimated survey-specific variation (random effect)
}

transformed parameters {
}

model {
  vector[K-1] cpr;       // cumulative probabilties
  vector[K]   pr;        // probabilities for each size class
  real tmp[J];           // probability of observed size class
  real mu;               // mean log fish size 
  real eps = 0.01;       // a small probability for random size class
  real c1;               // fraction of fish that do not fit distribution model
  real c2;               // probability randomly placed in size class

  c1 = (1.0 - eps);
  c2 = eps/K;
  
  // beta_0   ~ normal(3,0.5); // location-specific prior size (log)
  // beta_yr  ~ normal(0,0.1); // location-specific prior slope
  // beta_mpa ~ normal(0,0.1); // prior mpa effect
  yr_RE    ~ normal(0.0, sigma_yr);  // random interannual differences
  srv_RE   ~ normal(0.0, sigma_srv); // random survey differences

  for (i in 1:N) { // for each survey 
    // calculate mean body size for the survey, given location, year, and survey
    mu = beta_0[i_loc[i]] + beta_mpa*mpa[i] + beta_yr[i_loc[i]]*x[i] + 
         yr_RE[i_yr[i]] + srv_RE[i];

    // calculate cumulative probabilities of observing fish in each size class
    for (k in 1:(K-1)) {
      cpr[k] = normal_cdf(cutoff[k], mu, sigma_size); // cumulative probability
    }
    // calculate probabilities of observing fish in each size class (with random eps)
    pr[1] = c1*cpr[1] + c2; // probability observe in smallest size class
    for (k in 2:(K-1)) {
      pr[k] = c1*(cpr[k] - cpr[k-1]) + c2; // intermediate size classes
    }
    pr[K] = c1*(1.0 - cpr[K-1]) + c2; // probability of observing in largest size class
    
    for (k in 1:K) {
      target += y[i,k]*log(pr[k]); // add the log-likelihood term for each size class
    }
  }
}

generated quantities {
// int<lower=0, upper=1> y_new[T_new];
// for (t in 1:T_new){
//   y_new[t] = bernoulli_rng(p);
// }
}
