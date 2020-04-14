
/* 
   From https://arxiv.org/pdf/1707.01694.pdf

   From rstanarm
{ // hs+
    real log_half = -0.693147180559945286;
    target += normal_lpdf(z_beta | 0, 1);
    target += normal_lpdf(local[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
    target += normal_lpdf(local[3] | 0, 1) - log_half;
    // unorthodox useage of prior_scale as another df hyperparameter
    target += inv_gamma_lpdf(local[4] | 0.5 * prior_scale, 0.5 * prior_scale);
    target += normal_lpdf(global[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
    target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
  }
*/

data {
  int<lower=0> n; // number of observations
  int<lower=0> d; // number of predictors
  vector[n] y; // outputs
  matrix[n,d] x; // inputs
  real<lower=0> scale_icept; // prior std for the intercept
  real<lower=0> scale_global; // scale for the half -t prior for tau
  real<lower=1> nu_global; // degrees of freedom for the half -t priors for tau
  real<lower=1> nu_local; // degrees of freedom for the half - t priors for lambdas
  real<lower=0> slab_scale; // slab scale for the regularized horseshoe
  real<lower=0> slab_df; // slab degrees of freedom for the regularized horseshoe
}
parameters {
  real logsigma;
  real beta0;
  vector[d] z;
  real<lower=0> tau; // global shrinkage parameter
  vector<lower=0>[d] lambda; // local shrinkage parameter
  real<lower=0> caux;
}

transformed parameters {
  real<lower=0> sigma; // noise std
  vector<lower=0>[d] lambda_tilde; // ’ truncated ’ local shrinkage parameter
  real<lower=0> c; // slab scale
  vector[d] beta; // regression coefficients
  vector[n] f; // latent function values
  sigma = exp(logsigma);
  c = slab_scale * sqrt(caux);
  lambda_tilde = sqrt(c ^ 2 * square(lambda ) ./ (c ^ 2 + tau ^ 2 * square(lambda )) );
  beta = z .* lambda_tilde * tau;
  f = beta0 + x * beta;
}
model {
  // half -t priors for lambdas and tau , and inverse - gamma for c ^ 2
  z ~ normal(0, 1);
  lambda ~ student_t(nu_local, 0, 1);
  tau ~ student_t(nu_global, 0, scale_global * sigma);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  beta0 ~ normal(0, scale_icept);
  y ~ normal(f, sigma);
}
