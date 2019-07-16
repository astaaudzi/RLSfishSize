// Shane A Richards 6/11/2018
// includes a random effect associated with location, year, and survey
// combines identical observations

data {
  int  <lower = 1> N; // number of surveys
  int  <lower = 1> J; // number of locations
  int  <lower = 1> K; // number of size classes
  int  <lower = 1> L; // number of year classes
  real <lower = 0>            cutoff[K-1]; // size class cut-offs
  matrix[N,K]                 y;           // observed fish size classes
  real                        x[N];        // observed predictor variable
  int  <lower = 1, upper = J> gloc[N];     // observed location
  int  <lower = 1, upper = L> yr[N];       // observed year (index only)
}

parameters {
  real <lower =  2.0,   upper = 5.0> beta0;       // mean (log)size
  real <lower = -0.50,  upper = 0.5> beta1;       // predictor of mean (log)size
  real <lower =  0.10,  upper = 0.30> sigma_size;  // variation in (log)sizes
  real <lower =  0.001, upper = 1.0> sigma_gloc;  // variation among locations
  real <lower =  0.001, upper = 1.5> sigma_yr;    // variation among years
  real <lower =  0.001, upper = 1.0> sigma_srv;   // variation among years
  real gloc_RE[J];  // estimated location-specific variation
  real yr_RE[L];    // estimated year-specific variation
  real srv_RE[N];   // estimated survey-specific variation
}

transformed parameters {
}

model {
  vector[K-1] cpr;
  vector[K] pr;
  // real cpr[K-1] = 0.0;   // cumulative normal probability
  // real pr[K] = 0.0;      // probability of observed size class
  real mu;         // mean log fish size 
  real eps = 0.01; // a small probability
  real c1;
  real c2;

  c1 = (1.0 - eps);
  c2 = eps/K;
  
  // generate location and year specific variation 
  gloc_RE ~ normal(0.0, sigma_gloc); 
  yr_RE   ~ normal(0.0, sigma_yr); 
  srv_RE  ~ normal(0.0, sigma_srv); 

  for (i in 1:N) { // for each survey 
    mu = beta0 + beta1*x[i] + gloc_RE[gloc[i]] + yr_RE[yr[i]] + srv_RE[i];
    
    for (k in 1:(K-1)) {
      cpr[k] = normal_cdf(cutoff[k], mu, sigma_size);
    }
    
    pr[1] = c1*cpr[1] + c2;
    for (k in 2:(K-1)) {
      pr[k] = c1*(cpr[k] - cpr[k-1]) + c2;
    }
    pr[K] = c1*(1.0 - cpr[K-1]) + c2;
    
    for (k in 1:K) {
      target += y[i,k]*log(pr[k]); // add the log-likelihood term
    }
  }
}
