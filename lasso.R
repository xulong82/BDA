# l1/2 as loss functions

# l1 or least absolute deviation minimizes the sum of the absolute differences
# l2 or least squares error minimizs the sum of the square of the differences

# l1/2 as regularization

# regularization is a very important technique in machine learning to prevent overfitting
# difference between the l1 and l2 is that l2 is the sum of the square of the weights, while L1 is just the sum of the weights

library(glmnet)

rm(list = ls())
load("~/Library/R/3.4/library/glmnet/data/QuickStartExample.RData")

fit = glmnet(x, y, family = "gaussian")
coef(fit)

dim(x)

str(fit)

cvfit = cv.glmnet(x, y, family = "gaussian")
coef(cvfit)
