data {
  int<lower=1> g;  // number of genes
  int<lower=0> n;  // number of features
  int<lower=0> nnz; //number of non-zero entries in the feature matrix
  vector[g] B;  // Bayes Factors
  vector[nnz] w; //non-zero entries in the feature matrix
  int v[nnz];//column indices
  int u[g+1];//row start indices
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
  vector[g] pvec = inv_logit(csr_matrix_times_vector(g,n,w,v,u,Beta)+Beta0);
  Beta0 ~ normal(-1,4);
  Beta ~ normal(0,4);
  target+=sum(log( pvec .* B + (1 - pvec)));
}