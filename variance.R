# there is no error terms in logistic regression models, because there is no error terms in binomial distribution model
# but there are still residuals, variances, and such on

# Frank Harrell
# There is no error term in the Bernoulli distribution, there's just an unknown probability. The logistic model is a probability model.

# gaussian model
x = rnorm(100, 0, 1)
y = rnorm(100, 0, 1)

z1 = lm(y ~x)
z2 = glm(y ~x)

summary(z1)
sd(z1$residuals)
sd(z1$residuals)^2
var(z1$residuals)
var(z2$residuals)

summary(z2)
sum(z2$residuals^2)
z2$deviance
var(z2$residuals)

z2$deviance / (100 - 1)

# logistic model
x = rbinom(n = 100, size = 1, prob = 0.8)
z2 = glm(x ~ 1, family = "binomial")
summary(z2)

z2$fitted.values
z2$linear.predictors

log(z2$fitted.values / (1 - z2$fitted.values))

# deviance residual
residuals(z2)
sqrt(-2 * log(0.82))
-sqrt(-2 * log(1-0.82))

# partial residual
z2$residuals
(1 - 0.82) / 0.82 / 0.18
(0 - 0.82) / 0.82 / 0.18
