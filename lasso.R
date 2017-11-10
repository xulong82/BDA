# how to understand the parameter path?
# what are L1 and L2 regularizations for lasso and ridge?

library(glmnet)

rm(list = ls())
load("~/Library/R/3.4/library/glmnet/data/QuickStartExample.RData")

fit = glmnet(x, y, family = "gaussian")
coef(fit)

dim(x)

str(fit)

cvfit = cv.glmnet(x, y, family = "gaussian")
coef(cvfit)
