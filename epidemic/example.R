# propensity score matching is a useful tool for reducing selection bias and strengthening causal conclusions

library(dplyr)
library(ggplot2)
library(MatchIt)

mydata <- read.csv ("~/GitHub/stats/epidemic/newyork.csv") 
attach(mydata)
head(mydata)
str(mydata)

table(stw)

match = matchit(stw ~ tot + min + dis, data = mydata, method = "nearest", ratio = 1)
summary(match)
plot(match, type = "jitter")
plot(match, type = "hist")

matchdata = match.data(match)

glm(stw ~ tot, data = mydata, family = "binomial") %>% summary
glm(stw ~ tot, data = matchdata, family = "binomial") %>% summary

glm(stw ~ tot, data = mydata, family = "binomial") %>% summary

# propensity score
# use a logistic regression model to estimate the probability that a person received treatment
# this probability is the propensity score

fit = glm(stw ~ tot + min + dis, data = mydata, family = "binomial")
pscore = predict(fit, type = "response")

mydata$pscore = pscore
mydata$school = as.character(mydata$school)

# use propensity score to match samples

treat = mydata %>% filter(stw == 1)
control = mydata %>% filter(stw == 0)

school = NA

for(i in 1:nrow(treat)) {
  pscore.t = treat$pscore[i]
  index1 = which.min(abs(pscore.t - control$pscore)) # the nearest approach
  school1 = control$school[index1]
  control = control[-index1, ]
  school = c(school, school1)
}

school = school[-1]

mymatch = mydata[match(school, mydata$school), ]

mymatch = rbind(mymatch, treat)

setdiff(mymatch$school, as.character(matchdata$school))
setdiff(as.character(matchdata$school), mymatch$school)

plot(mymatch$pscore, jitter(mymatch$stw), xlim = c(0, 0.4), ylim = c(-1, 3))

# note that in the above matching method, we dropped many samples in the control group

# there are other ways to use propensity score that avoid this drawback

# 1. use propensity score as covariate, so we use all samples to estimate tx effects
# 2. use propensity score as weights to the response, again we use all samples

# outcome analysis w/o adjustment
lm(re78 ~ treat + black + hispan + married, data = lalonde) %>% summary

# balance analysis

lm(re74 ~ treat, data = lalonde) %>% summary
glm(nodegree ~ treat, data = lalonde, family = "binomial") %>% summary

# estimate propensity score

pscore = glm(treat ~ age + educ + nodegree + re74 + re75, data = lalonde, family = binomial)
summary(pscore)

lalonde$pscore = predict(pscore, type = "response")

# weight estimation using propensity score

lalonde$weight = ifelse(lalonde$treat == 1, 1/lalonde$pscore, 1/(1-lalonde$pscore))

lm(re74 ~ treat, data = lalonde, weights = weight) %>% summary

# use propensity score as weight in outcome analysis

lm(re78 ~ treat + black + hispan + married, data = lalonde) %>% summary
lm(re78 ~ treat + black + hispan + married, data = lalonde, weights = weight) %>% summary

# in a weighted fit, less weight is given to the less precise measurements 
# and more weight to more precise measurements when estimating the unknown parameters in the model

lalonde$re78.adj = lalonde$re78 * lalonde$weight
lm(re78.adj ~ treat + black + hispan + married, data = lalonde) %>% summary

# it isn't simply multipying weight to the reponse, instead implemented in the object function of parameter estimation process

# \sum wi*(yi - \hat{yi})^2

# weighting method may lose power in cases where 
# the variables under considerations simply reflect selection bias, but not confouding in estimating treatment effects
