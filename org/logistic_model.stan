data {
  int<lower=1> g;  // number of data points
  int<lower=0> n;  // number of clusters
  vector[g] B;  // Bayes Factors
  matrix[g,n] A; //Annotation Matrix
}
parameters {
  vector[n] Beta; // Effect sizes
  real Beta0; //intercept
}
/* transformed parameters { */
/* vector<lower=0,upper=1>[G] pvec=inv_logit(A * Beta + Beta0); */
/* } */
model {
  // likelihood
  vector[g] pvec = inv_logit((A*Beta)+Beta0);
  Beta0 ~ normal(-1,4);
  Beta ~ normal(0,4);
  target+=sum(log( pvec .* B + (1 - pvec)));
}