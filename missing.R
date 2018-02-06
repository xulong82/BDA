# there are 3 usual practices to handle missing values in predictors under the linear regression framework: 
# 1. delete the samples; this is the last choice because we lose information of predictors that do not have missing fields
# 2. imputation; this is a reasonable approach for ordinal predictors, under the missing at random assumption
# 3. categorize the variables: this is the most desirable approach for categorical predictor

# what we will get by treating missing fields of categorical predictors as a separate category?

# with simulations, we found by introducing a new "missing" category:
# 1. mean estimands of known categories are not affected 
# 2. variance estimands of known categories are affected

# how variance estimands are affected?
# turns out lm computes a combind variances for all categories of a categorical variables
# instead of separate variance estimands for each category
# this combined variance estimand was further weighted by sample size of each category to derive estimands of sd.error
# with this, if the missing fields are missed at random, the combined variance estimand will not be affected
# sd.error and p-value estimands of each known categories are also not affected
# if the missing fields have larger variance (e.g. mostly come from a group of large variance) 
# the combined variance estimand will be larger, and consequently st.error and p-values of known categories will be larger
# vice versa

library(dplyr)

x = c(rep("a", 1e2), rep("b", 1e2), rep("c", 1e2))
x = as.factor(x)

y = c(rnorm(1e2, 1, 1), rnorm(1e2, 2, 2), rnorm(1e2, 3, 3))

lm(y ~ x - 1) %>% summary

x.miss = rep("m", 1e1)

y.miss = sample(y, 1e1) # missing at random

y.miss = sample(y[x == "a"], 1e1) # skewed to one level
y.miss = sample(y[x == "b"], 1e1) # skewed to one level
y.miss = sample(y[x == "c"], 1e1) # skewed to one level

var(y[x == "a"])
var(y[x == "b"])
var(y[x == "c"])

var(y.miss)

x = c(rep("a", 1e2), rep("b", 1e2), rep("c", 1e2))
x.new = factor(c(x, x.miss))

y.new = c(y, y.miss)

lm(y.new ~ x.new - 1) %>% summary

# end

# binary indicator?

x1 <- rnorm(1e2)
y1 <- 2*x1 + rnorm(1e2, 3, 1)

x2 <- rnorm(30)
y2 <- rep(0, 30)

yy <- c(y1, y2)
xx <- c(x1, x2)

i <- c(rep(1, 1e2), rep(0, 30))

zz = i * xx

hist(xx, n = 1e2)
hist(yy, n = 1e2)

plot(zz, yy)

lm(yy ~ 0 + i + zz) %>% summary
lm(yy ~ 0 + i + i:xx) %>% summary
lm(yy[1:1e2] ~ 0 + i[1:1e2] + xx[1:1e2]) %>% summary

# binary indicator approach does not affect mean estimands, however, tends to shrink the std.err estimands
