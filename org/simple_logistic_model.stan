
  data {
    int<lower=1> G;  // number of data points
    int<lower=0> F;  // number of clusters
    vector[G] B;  // Bayes Factors
    matrix[G,F] A; //Annotation Matrix
  }
  parameters {
    vector[F] Beta; // Effect sizes
    real Beta0; //intercept
  }
  /* transformed parameters { */
  /* vector<lower=0,upper=1>[G] pvec=inv_logit(A * Beta + Beta0); */
  /* } */
  model {
    // likelihood
    Beta0 ~ normal(-3.85,5);
for
    Beta ~ normal(0,4.5);
    for(n in 1:G){
      real pv=Beta0;
      for(k in 1:F){
        pv+=A[n,k]*Beta[k];
      }
      target+=log( inv_logit(pv) * B[n] + (1 - inv_logit(pv)));
    }
  }
