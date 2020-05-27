
data {
  int<lower=1> G;  // number of data points
  int<lower=0> F;  // number of clusters
  vector[G] B;  // Bayes Factors
  matrix[G,F] A; //Annotation Matrix
}
parameters {
  vector[F] Beta; // Effect sizes
  real Beta0; //intercept
  real<lower=0> shape;
  real<lower=0> rate;
}
/* transformed parameters { */
/* vector<lower=0,upper=1>[G] pvec=inv_logit(A * Beta + Beta0); */
/* } */
model {
  // likelihood
  vector[G] pvec = inv_logit(A*Beta+Beta0);
  Beta0 ~ normal(-2,5);
  Beta ~ gamma(shape,rate);
  target+=sum(log( pvec .* B + (1 - pvec)));
}
