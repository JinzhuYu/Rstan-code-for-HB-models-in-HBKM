
## RStan code for the HB models used in HBKM and HBRM

  # HBKM in stan

  Poisson_HBK='
  data {
  int<lower=1> I; 
  int<lower=1> J;  
  int<lower=1> Y[I,J]; 
  }
  
  parameters {
  real<lower=0> alpha[J];
  real<lower=0> beta[J];
  matrix<lower=0>[I,J] lambda; 
  }
  
  model {
  for(i in 1:I){
  for (j in 1:J){
  target +=cauchy_lpdf(alpha[j]|0,5); 
  target +=cauchy_lpdf(beta[j]|0,5);
  target +=gamma_lpdf(lambda[i,j]|alpha[j],beta[j]);
  target +=poisson_lpmf(Y[i,j]|lambda[i,j]);
  }
  }
  }
  
  '

  # HBRM in stan
  Poisson_HBR='
  data {
  int<lower=1> I; 
  int<lower=1> J;  
  int<lower=1> y[I,J];   
  real<lower=1> x[I,J];  
  }
  
  parameters {
  vector[2] mu;
  matrix[2,J] theta;
  corr_matrix[2] Omega;  
  vector<lower=0>[2] tau;  
  }
  
  model {
  matrix[I,J] lambda;
  mu~normal(0,5);
  tau ~ cauchy(0,5);  
  Omega ~ lkj_corr(2);   
  for(i in 1:I){
  for (j in 1:J){
  theta[,j]~multi_normal(mu, quad_form_diag(Omega, tau));
  lambda[i,j] = exp(theta[1,j]*x[i,j]+theta[2,j]);
  y[i,j]~poisson(lambda[i,j]);
  }
  }
  }
  
  '

library(rstan)

n_sam=2000
n_warmup=n_sam/2
n_chain=4

# fit HBKM
Poisson_HBK_dat=list(Y=Y_train,I=nrow(Y_train),J=jj)  # J=5, which is the number of groups in the dataset
fit_HBK=stan(model_code=Poisson_HBK, data = Poisson_HBK_dat,iter =(n_sam+n_warmup), warmup=n_warmup,chains = n_chain)
fitresult_HBK= extract(fit_HBK, permuted = TRUE)

# fit HBRM
Poisson_HBR_dat=list(y=rbind(Y_train,Y_tune),x=rbind(X_train,X_tune),I=length(c(Y_train[,1],Y_tune[,1])),J=jj)
fit_HBR=stan(model_code=Poisson_HBR, data = Poisson_HBR_dat,iter =(n_sam+n_warmup), warmup=n_warmup,chains = n_chain)
fitresult_HBR= extract(fit_HBR, permuted = TRUE)
