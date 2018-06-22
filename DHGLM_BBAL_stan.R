###### pacakges ######
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###### data ######
data <- read.table("bugsdata5.txt")
no_of_cores <- 3
### removing NAs
I <-is.na(data$ids)|
  is.na(data$Sexe)|
  is.na(data$Cycle)|
  is.na(data$CycleBand)|
  is.na(data$Rsnum)|
  is.na(data$RSplus1num)|
##  is.na(data$CyclePairBand)|
  is.na(data$Stade)|
  is.na(data$ntrajet)|
  is.na(data$Tripdur)|
  is.na(data$Age)|
  is.na(data$IODAnnual)| 
  is.na(data$SOIAnnual)|  
  is.na(data$AAOAnnual) ##|  
##  is.na(data$ONIAnnual)
 
data <- data[!I,]

### scaling and centering
## data$ntrajet<- as.numeric(scale(data$ntrajet))
data$Tripdur <- as.numeric(scale(sqrt(data$Tripdur)))
data$Age<- as.numeric(scale(data$Age))
data$IODAnnual<- as.numeric(scale(data$IODAnnual))
data$SOIAnnual<- as.numeric(scale(data$SOIAnnual))
data$AAOAnnual<- as.numeric(scale(data$AAOAnnual))
##data$ONIAnnual<- as.numeric(scale(data$ONIAnnual))

db <- data.frame(Cycle=sort(unique(data$Cycle)),Cycle.n=1:length(unique(data$Cycle)))
data$Cycle.n <- db$Cycle.n[match(data$Cycle,db$Cycle)]

Xm <- unname( model.matrix(~ 1 + Age + (Sexe + Stade) * IODAnnual, data))
Xv <- unname( model.matrix(~ 1 + Age + (Sexe + Stade) * IODAnnual, data))

##### creating BUGS datalist #####
stan.data <- list(
    Km = ncol(Xm),
    Kv = ncol(Xv),
    N = nrow(Xm),
    nid = length(unique(data$ids)),
    npr = length(unique(data$PairID)),
    nyr = length(unique(data$Cycle)),
    id = as.numeric(droplevels(data$ids)),
    pr = as.numeric(droplevels(as.factor(data$PairID))),
    yr = data$Cycle.n,
    Xm = Xm,
    Xv = Xv,
    Y = data$Tripdur
)


sink("Stan_Model.stan")
cat("
data {
  int<lower=1>    Km; // number of fixed effects on mean
  int<lower=1>    Kv; // number of fixed effects on variance
  int<lower=0>    N; // number of observations
  int<lower=1>    nid; // number of individuals
  int<lower=1>    npr;
  int<lower=1>    nyr;
  int<lower=1, upper=nid> id[N]; // individual identity
  int<lower=1, upper=npr> pr[N]; // group
  int<lower=1, upper=nyr> yr[N];
  vector[Km]   Xm[N]; // Fixed effects design matrix on mean
  vector[Kv]   Xv[N];
  real         Y[N]; // response variable

}
parameters {
  row_vector[Km] beta; // fixed effects
 row_vector[Kv] gamma; // fixed effects

  vector<lower=0>[2] sigma_id;
  vector<lower=0>[2] sigma_pr;
  vector<lower=0>[2] sigma_yr;

  cholesky_factor_corr[2] L_id;
  cholesky_factor_corr[2] L_pr;
  cholesky_factor_corr[2] L_yr;

  matrix[2,nid] z_id;
  matrix[2,npr] z_pr;
  matrix[2,nyr] z_yr;

}
transformed parameters{
   matrix[2, nid] uid = diag_pre_multiply(sigma_id, L_id) * z_id; 
   matrix[2, npr] upr = diag_pre_multiply(sigma_pr, L_pr) * z_pr;
   matrix[2, nyr] uyr = diag_pre_multiply(sigma_yr, L_yr) * z_yr; 
}
model {
    vector[N] mu;
    vector[N] sigma_e;

    for(n in 1:N) {
    mu[n] = beta * Xm[n] + uid[1, id[n]] + upr[1,pr[n]] + + uyr[1, yr[n]];
    sigma_e[n] = exp(gamma * Xv[n] + uid[2, id[n]] + upr[2,pr[n]] + uyr[2, yr[n]]);

}

      Y ~ normal(mu, sigma_e);
//priors
    to_vector(beta) ~ normal(0, 1);
    to_vector(gamma) ~ normal(0, 1);
    to_vector(z_id) ~ normal(0, 1);
    to_vector(z_pr) ~ normal(0, 1);
    to_vector(z_yr) ~ normal(0, 1);
    L_id ~ lkj_corr_cholesky(1);
    sigma_id ~ cauchy(0, 5);
    L_pr ~ lkj_corr_cholesky(1);
    sigma_pr ~ cauchy(0, 5);
    L_yr ~ lkj_corr_cholesky(1);
    sigma_yr ~ cauchy(0, 5);

}
generated quantities {
    cov_matrix[2] Sigma_id;
    cov_matrix[2] Sigma_pr;
    cov_matrix[2] Sigma_yr;
    corr_matrix[2] cor_id;
    corr_matrix[2] cor_pr;
    corr_matrix[2] cor_yr;
    real R_int; 

    Sigma_id = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_id, L_id));
    Sigma_pr = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_pr, L_pr));
    Sigma_yr = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_yr, L_yr));

    cor_id = multiply_lower_tri_self_transpose(L_id);
    cor_pr = multiply_lower_tri_self_transpose(L_pr);
    cor_yr = multiply_lower_tri_self_transpose(L_yr);

    R_int = sigma_id[1]^2 / (sigma_id[1]^2 + sigma_pr[1]^2 + sigma_yr[1]^2 + gamma[1]^2);

}
    ",fill = TRUE)
  sink()


dh_nog<- stan(file="Stan_Model.stan", data = stan.data, iter = 4000, warmup = 2000, thin = 1, chains = 4) #, control = list(adapt_delta = 0.95))
##save(dh_nog,file="dhglm_nog.rda")
