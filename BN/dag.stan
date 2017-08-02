data {
  
  int Ntrain;
  int Ntest;
  
  int<lower=0, upper=1> sVirgi[Ntrain];
  int<lower=0, upper=1> sVersi[Ntrain];
  
  vector[Ntrain] sLength1;
  vector[Ntrain] sWidth1;
  vector[Ntrain] pLength1;
  vector[Ntrain] pWidth1;
  
  
  vector[Ntest] sLength;
  vector[Ntest] sWidth;
  vector[Ntest] pLength;
  vector[Ntest] pWidth;
  
}

transformed data {
  vector[Ntrain] sVirgiVec;
  vector[Ntrain] sVersiVec;
  vector[Ntrain] on;
  vector[Ntrain] off;
  
  for (i in 1:Ntrain) {
    sVirgiVec[i] <- sVirgi[i] - 0.0;
    sVersiVec[i] <- sVersi[i] - 0.0;
    on[i] <- 1.0;
    off[i] <- 0.0;
  }
  
}

parameters {
  
  real sVir_Intercept;
  real sVer_Intercept;
  
  real pWidth_Intercept;
  real pWidth_sVir;
  real pWidth_sVer;
  real<lower=0> pWidth_sigma;
  
  real sLength_Intercept;
  real sLength_pWidth;
  real<lower=0> sLength_sigma;
  
  real sWidth_Intercept;
  real sWidth_sVir;
  real sWidth_sVer;
  real sWidth_pWidth;
  real sWidth_sLength;
  real<lower=0> sWidth_sigma;
  
  real pLength_Intercept;
  real pLength_sVir;
  real pLength_sVer;
  real pLength_pWidth;
  real pLength_sLength;
  real<lower=0> pLength_sigma;
  
}

model {
  
  sVersi ~ bernoulli_logit(sVer_Intercept);
  sVirgi ~ bernoulli_logit(sVir_Intercept);
  sVer_Intercept ~ normal(0, 100);
  sVir_Intercept ~ normal(0, 100);
  
  pWidth1 ~ normal(pWidth_Intercept + pWidth_sVer*sVersiVec + pWidth_sVir*sVirgiVec, pWidth_sigma);
  pWidth_Intercept ~ normal(0, 100);
  pWidth_sVer ~ normal(0, 100); 
  pWidth_sVir ~ normal(0, 100);
  pWidth_sigma ~ cauchy(0 ,5);
  
  sLength1 ~ normal(sLength_Intercept + sLength_pWidth*pWidth1, sLength_sigma);
  sLength_Intercept ~ normal(0, 100);
  sLength_pWidth ~ normal(0, 100);
  sLength_sigma ~ cauchy(0, 5);
  
  sWidth1 ~ normal(sWidth_Intercept + sWidth_sVer*sVersiVec + sWidth_sVir*sVirgiVec +
                  sWidth_pWidth*pWidth1 + sWidth_sLength*sLength1, sWidth_sigma);
  sWidth_Intercept ~ normal(0, 100);
  sWidth_sVer ~ normal(0, 100); 
  sWidth_sVir ~ normal(0, 100);
  sWidth_sigma ~ cauchy(0 ,5);
  
  pLength1 ~ normal(pLength_Intercept + pLength_sVer*sVersiVec + pLength_sVir*sVirgiVec +
                  pLength_pWidth*pWidth1 + pLength_sLength*sLength1, pLength_sigma);
  pLength_Intercept ~ normal(0, 100);
  pLength_sVer ~ normal(0, 100); 
  pLength_sVir ~ normal(0, 100);
  pLength_sigma ~ cauchy(0 ,5);
  
}

generated quantities {
  vector[Ntest] virgi;
  vector[Ntest] notVirgi;
  vector[Ntest] versi;
  vector[Ntest] notVersi;
  
  for (i in 1:Ntest){ 
    
  virgi[i] <- exp(log(inv_logit(sVir_Intercept)) +
           log(1-inv_logit(sVer_Intercept)) +
           normal_log(pWidth[i], pWidth_Intercept + pWidth_sVir*on[i] + pWidth_sVer*off[i], pWidth_sigma) +
           normal_log(sLength[i], sLength_Intercept + pWidth[i]*sLength_pWidth, sLength_sigma) +
           normal_log(sWidth[i], sWidth_Intercept  + sWidth_sVir*on[i] + sWidth_sVer*off[i] + 
                  sLength[i]*sWidth_sLength + pWidth[i]*sWidth_pWidth, sWidth_sigma) +
           normal_log(pLength[i], pLength_Intercept + pLength_sVir*on[i] + pLength_sVer*off[i] +
                  pWidth[i]*pLength_pWidth + sLength[i]*pLength_sLength, pLength_sigma)) + 
           exp(log(inv_logit(sVir_Intercept)) +
           log(inv_logit(sVer_Intercept)) +
           normal_log(pWidth[i], pWidth_Intercept + pWidth_sVir*on[i] + pWidth_sVer*on[i], pWidth_sigma) +
           normal_log(sLength[i], sLength_Intercept + pWidth[i]*sLength_pWidth, sLength_sigma) +
           normal_log(sWidth[i], sWidth_Intercept  + sWidth_sVir*on[i] + sWidth_sVer*on[i] + 
                  sLength[i]*sWidth_sLength + pWidth[i]*sWidth_pWidth, sWidth_sigma) +
           normal_log(pLength[i], pLength_Intercept + pLength_sVir*on[i] + pLength_sVer*on[i] +
                  pWidth[i]*pLength_pWidth + sLength[i]*pLength_sLength, pLength_sigma));
                  
    notVirgi[i] <- exp(log(1-inv_logit(sVir_Intercept)) +
           log(1-inv_logit(sVer_Intercept)) +
           normal_log(pWidth[i], pWidth_Intercept + pWidth_sVir*off[i] + pWidth_sVer*off[i], pWidth_sigma) +
           normal_log(sLength[i], sLength_Intercept + pWidth[i]*sLength_pWidth, sLength_sigma) +
           normal_log(sWidth[i], sWidth_Intercept  + sWidth_sVir*off[i] + sWidth_sVer*off[i] + 
                  sLength[i]*sWidth_sLength + pWidth[i]*sWidth_pWidth, sWidth_sigma) +
           normal_log(pLength[i], pLength_Intercept + pLength_sVir*off[i] + pLength_sVer*off[i] +
                  pWidth[i]*pLength_pWidth + sLength[i]*pLength_sLength, pLength_sigma)) + 
           exp(log(1-inv_logit(sVir_Intercept)) +
           log(inv_logit(sVer_Intercept)) +
           normal_log(pWidth[i], pWidth_Intercept + pWidth_sVir*off[i] + pWidth_sVer*on[i], pWidth_sigma) +
           normal_log(sLength[i], sLength_Intercept + pWidth[i]*sLength_pWidth, sLength_sigma) +
           normal_log(sWidth[i], sWidth_Intercept  + sWidth_sVir*off[i] + sWidth_sVer*on[i] + 
                  sLength[i]*sWidth_sLength + pWidth[i]*sWidth_pWidth, sWidth_sigma) +
           normal_log(pLength[i], pLength_Intercept + pLength_sVir*off[i] + pLength_sVer*on[i] +
                  pWidth[i]*pLength_pWidth + sLength[i]*pLength_sLength, pLength_sigma));
                  
    versi[i] <- exp(log(inv_logit(sVir_Intercept)) +
           log(inv_logit(sVer_Intercept)) +
           normal_log(pWidth[i], pWidth_Intercept + pWidth_sVir*on[i] + pWidth_sVer*on[i], pWidth_sigma) +
           normal_log(sLength[i], sLength_Intercept + pWidth[i]*sLength_pWidth, sLength_sigma) +
           normal_log(sWidth[i], sWidth_Intercept  + sWidth_sVir*on[i] + sWidth_sVer*on[i] + 
                  sLength[i]*sWidth_sLength + pWidth[i]*sWidth_pWidth, sWidth_sigma) +
           normal_log(pLength[i], pLength_Intercept + pLength_sVir*on[i] + pLength_sVer*on[i] +
                  pWidth[i]*pLength_pWidth + sLength[i]*pLength_sLength, pLength_sigma)) + 
           exp(log(1-inv_logit(sVir_Intercept)) +
           log(inv_logit(sVer_Intercept)) +
           normal_log(pWidth[i], pWidth_Intercept + pWidth_sVir*off[i] + pWidth_sVer*on[i], pWidth_sigma) +
           normal_log(sLength[i], sLength_Intercept + pWidth[i]*sLength_pWidth, sLength_sigma) +
           normal_log(sWidth[i], sWidth_Intercept  + sWidth_sVir*off[i] + sWidth_sVer*on[i] + 
                  sLength[i]*sWidth_sLength + pWidth[i]*sWidth_pWidth, sWidth_sigma) +
           normal_log(pLength[i], pLength_Intercept + pLength_sVir*off[i] + pLength_sVer*on[i] +
                  pWidth[i]*pLength_pWidth + sLength[i]*pLength_sLength, pLength_sigma));
                  
    notVersi[i] <- exp(log(1-inv_logit(sVir_Intercept)) +
           log(1-inv_logit(sVer_Intercept)) +
           normal_log(pWidth[i], pWidth_Intercept + pWidth_sVir*off[i] + pWidth_sVer*off[i], pWidth_sigma) +
           normal_log(sLength[i], sLength_Intercept + pWidth[i]*sLength_pWidth, sLength_sigma) +
           normal_log(sWidth[i], sWidth_Intercept  + sWidth_sVir*off[i] + sWidth_sVer*off[i] + 
                  sLength[i]*sWidth_sLength + pWidth[i]*sWidth_pWidth, sWidth_sigma) +
           normal_log(pLength[i], pLength_Intercept + pLength_sVir*off[i] + pLength_sVer*off[i] +
                  pWidth[i]*pLength_pWidth + sLength[i]*pLength_sLength, pLength_sigma)) + 
           exp(log(inv_logit(sVir_Intercept)) +
           log(1-inv_logit(sVer_Intercept)) +
           normal_log(pWidth[i], pWidth_Intercept + pWidth_sVir*on[i] + pWidth_sVer*off[i], pWidth_sigma) +
           normal_log(sLength[i], sLength_Intercept + pWidth[i]*sLength_pWidth, sLength_sigma) +
           normal_log(sWidth[i], sWidth_Intercept  + sWidth_sVir*on[i] + sWidth_sVer*off[i] + 
                  sLength[i]*sWidth_sLength + pWidth[i]*sWidth_pWidth, sWidth_sigma) +
           normal_log(pLength[i], pLength_Intercept + pLength_sVir*on[i] + pLength_sVer*off[i] +
                  pWidth[i]*pLength_pWidth + sLength[i]*pLength_sLength, pLength_sigma));

  }              

}