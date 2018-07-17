
data {
  int<lower=1>    Km; // number of fixed effects on mean
  int<lower=1>    Kv; // number of fixed effects on variance
  int<lower=0>    N; // number of observations
  int<lower=1>    nid; // number of individuals
  int<lower=1>    nyr;
  int<lower=1, upper=nid> id[N]; // individual identity
  int<lower=1, upper=nyr> yr[N];
  vector[Km]   Xm[N]; // Fixed effects design matrix on mean
  vector[Kv]   Xv[N];
  real         Y[N]; // response variable

}
parameters {
  row_vector[Km] beta; // fixed effects
  row_vector[Kv] gamma; // fixed effects

  vector<lower=0>[2] sigma_id;
  vector<lower=0>[2] sigma_yr;

  cholesky_factor_corr[2] L_id;
  cholesky_factor_corr[2] L_yr;

  matrix[2,nid] z_id;
  matrix[2,nyr] z_yr;

}
transformed parameters{
   matrix[2, nid] uid = diag_pre_multiply(sigma_id, L_id) * z_id; 
   matrix[2, nyr] uyr = diag_pre_multiply(sigma_yr, L_yr) * z_yr; 
}
model {
    vector[N] mu;
    vector[N] sigma_e;

    for(n in 1:N) {
    mu[n] = beta * Xm[n] + uid[1, id[n]] + uyr[1, yr[n]];
    sigma_e[n] = exp(gamma * Xv[n] + uid[2, id[n]] + uyr[2, yr[n]]);

}

      Y ~ normal(mu, sigma_e);
//priors
    to_vector(beta) ~ normal(0, 1);
    to_vector(gamma) ~ normal(0, 1);
    to_vector(z_id) ~ normal(0, 1);
    to_vector(z_yr) ~ normal(0, 1);
    L_id ~ lkj_corr_cholesky(1);
    sigma_id ~ cauchy(0, 5);
    L_yr ~ lkj_corr_cholesky(1);
    sigma_yr ~ cauchy(0, 5);

}
generated quantities {
    cov_matrix[2] Sigma_id;
    cov_matrix[2] Sigma_yr;
    corr_matrix[2] cor_id;
    corr_matrix[2] cor_yr;
    real R_int; 

    Sigma_id = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_id, L_id));
    Sigma_yr = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_yr, L_yr));

    cor_id = multiply_lower_tri_self_transpose(L_id);
    cor_yr = multiply_lower_tri_self_transpose(L_yr);

    R_int = sigma_id[1]^2 / (sigma_id[1]^2 + sigma_yr[1]^2 + gamma[1]^2);

}
    
