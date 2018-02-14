# l1/2 as loss functions

# l1 or least absolute deviation minimizes the sum of the absolute differences
# l2 or least squares error minimizs the sum of the square of the differences

# l1/2 as regularization

# regularization is a very important technique in machine learning to prevent overfitting
# difference between the l1 and l2 is that l2 is the sum of the square of the weights, while L1 is just the sum of the weights

# Rob Tribshrani lasso page: http://statweb.stanford.edu/~tibs/lasso/simple.html
# l1: minimize sum( (y-yhat)^2 ) subject to sum[absolute value(bj)] <= s

# efficient algorithms exists to solve the lasso solutions, such as least angle regression algorithm that follows:
# 1. start with all coefficients bj equal to zero
# 2. find the predictor xj most correlated with y
# 3. increase the coefficient bj in the direction of the sign of its correlation with y, 
# take residuals r=y-yhat along the way, stop when some other predictor xk has as much correlation with r as xj has

# issues with lasso
# 1. if you have correlated informative predictors, lasso tends to choose one and push the others to 0; 
# 2. lasso tends to choose more variables than a real model has; 
# 3. lasso selections are unstable, depending on the training sample; 
# 4. the betas are biased, and therefore cannot be used for a classical parametric hypothesis testing; 

# use lasso to make inferences on predictors
# bootstrap the data and count how many times each variable is selected, divide by number of resamples, 

library(glmnet)
library(dplyr)

rm(list = ls())
load("~/Library/R/3.4/library/glmnet/data/QuickStartExample.RData")

fit = glmnet(x, y, family = "gaussian")

plot(fit, xvar = "lambda", label = T)
plot(fit, xvar = "norm", label = T)
print(fit)

coef(fit)
coef(fit, s = 0.1)

cvfit = cv.glmnet(x, y, family = "gaussian")
coef(cvfit)

cvfit$lambda.min %>% log # gives minimum mean cross-validation error
cvfit$lambda.1se %>% log # gives the most regularized model such taht error is within one standard error of the minimum

coef(cvfit, s = "lambda.min")
coef(cvfit, s = "lambda.min") %>% as.matrix

# coef(cvfit) results are by taking the complete training data, is it?

test = glmnet(x, y, family = "gaussian", lambda = cvfit$lambda.min)
coef(test) %>% as.matrix

# ? why l1 penalty put coefficients of non-relevant predictors all the way to 0

# interaction terms

data <- data.frame(matrix(rnorm(9 * 10), ncol = 9))
names(data) <- c(paste0("x", 1:8), "y")

f <- as.formula(y ~ .*.) # all interactions
f <- as.formula(y ~ x1*.) # one with all others

y <- data$y

x <- model.matrix(f, data)[, -1]
glmnet(x, y)

