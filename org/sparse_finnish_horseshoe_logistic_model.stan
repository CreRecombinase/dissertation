  data {
    int<lower=1> g;  // number of genes
    int<lower=0> n;  // number of features
    int<lower=0> nnz; //number of non-zero entries in the feature matrix
    real m0 ;           // Expected number of large slopes
    vector[g] B;  // Bayes Factors
    vector[nnz] w; //non-zero entries in the feature matrix
    int v[nnz];//column indices
    int u[g+1];//row start indices

  }
  transformed data {
    real slab_scale = 3;    // Scale for large slopes
    real slab_scale2 = square(slab_scale);
    real slab_df = 25;      // Effective degrees of freedom for large slopes
    real half_slab_df = 0.5 * slab_df;
  }
  parameters {
    vector[n] beta_tilde; // Effect sizes
    real Beta0; //intercept
    vector<lower=0>[n] lambda;
    real<lower=0> c2_tilde;
    real<lower=0> tau_tilde;
    real<lower=0> sigma;
  }
transformed parameters{
  vector[n] beta;
  {
    real tau0 = (m0 / (n - m0)) * (sigma / sqrt(1.0 * g));
    real tau = tau0 * tau_tilde; // tau ~ cauchy(0, tau0)

    // c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scale2)
    // Implies that marginally beta ~ student_t(slab_df, 0, slab_scale)
    real c2 = slab_scale2 * c2_tilde;

    vector[n] lambda_tilde = sqrt( c2 * square(lambda) ./ (c2 + square(tau) * square(lambda)) );

    // beta ~ normal(0, tau * lambda_tilde)
    beta = tau * lambda_tilde .* beta_tilde;

  }

}
  model {

    // likelihood
    vector[g] pvec = inv_logit(csr_matrix_times_vector(g,n,w,v,u,beta)+Beta0);
    beta_tilde ~ normal(0, 1);
    Beta0 ~ normal(-1,2);
    lambda ~ cauchy(0, 1);
    tau_tilde ~ cauchy(0, 1);
    c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);
    target+=sum(log( pvec .* B + (1 - pvec)));
  }
