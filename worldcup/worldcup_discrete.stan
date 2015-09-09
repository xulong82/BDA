data {
  int nteams;
  int ngames;
  vector[nteams] prior_score;
  int team1[ngames];
  int team2[ngames];
  int score1[ngames];
  int score2[ngames];
  int maxrange;
  int nrange;
}
transformed data {
  vector[nrange] range;
  vector[ngames] dif;
  vector[ngames] sqrt_dif;
  for (i in 1:nrange)
    range[i] <- sqrt (abs (i - (maxrange + 1)));
  dif <- score1 - score2;
  for (i in 1:ngames)
    sqrt_dif[i] <- (step(dif[i]) - .5)*sqrt(fabs(dif[i]));
}
}
parameters {
  real b;
  real<lower=0> sigma_a;
  real<lower=0> sigma_y;
  vector[nteams] eta_a;
  real<lower=1> df;
}
transformed parameters {
  vector[nteams] a;
  a <- b*prior_score + sigma_a*eta_a;
}  
model {
  vector[nrange] probs;
  vector[nrange] probs_normalized;
  int dif_categorized;
  df ~ gamma(2,0.1);
  eta_a ~ normal(0,1);
  for (i in 1:ngames){
    for (j in 1:nrange){
      probs[j] <- exp(student_t_log(range[j], df, a[team1[i]]-a[team2[i]], sigma_y));
    }
    probs_normalized <- probs/sum(probs);
    dif_categorized <- score1[i] - score2[i] + (maxrange + 1);
    dif_categorized ~ categorical(probs_normalized);
  }
}
