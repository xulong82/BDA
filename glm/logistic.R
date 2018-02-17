# Where the prediction probabilities come from? 
# What do they mean?
# What does a prediction probability of 0.5 mean?
# Why the prediction probabilities in Atezo were overwhelmly low?
# How imbalanced response affects logistic regression?

library(dplyr)
rm(list = ls())

# logistic regression under the hood

# it model P(1) by linear combination of predictors via a link function
# P(1) = 1 / (1 + exp(-z))
# z = b0 + b1 * x

pfun = function(z) 1/(1+exp(-z))
z = rnorm(1e2, 0, 3)
plot(z, pfun(z))
abline(h = 0.5)
abline(v = 0.0)

# this is the sigmoid function, and it is nice in turning all real numbers to [0, 1] with a smooth curve
# with the sigmoid function, prediction probability of each subject is all up to the model parameters: b0, b1
# to estimate the model parameters is all about the cost function
# cost function of logistic model is the product of pfun(z) and (1-pfun(z)) for each subject
# we find b0 and b1 that maximize the cost function

# test using real-life data
library(mlbench)
data(BreastCancer, package="mlbench")

bc <- BreastCancer[complete.cases(BreastCancer), ] 
for(i in 2:10) bc[, i] <- as.numeric(as.character(bc[, i]))
bc$Class <- ifelse(bc$Class == "malignant", 1, 0)
bc$Class <- factor(bc$Class, levels = c(0, 1))

# let's stratify the subjects by case or control
bc0 = bc[bc$Class == 0, ] # all control
bc1 = bc[bc$Class == 1, ] # all cases

# what if all subjects were 0
fit <- glm(Class ~ Cell.shape, family = "binomial", data=bc0)
# the cost function does not converge, because the cost function becomes monotonical 
# and the predicted probabilities are all extremely low, which actually fit the data and the cost function
predict(fit, type = "response") %>% hist

# what if adding a few cases
data = rbind(bc0, bc1[sample(1:nrow(bc1), 3), ])
fit <- glm(Class ~ Cell.shape, family = "binomial", data=data)
summary(fit)
(fit.p = predict(fit, type = "response")) %>% hist
# the predicted probabilities are still very low for each subject
# note that final cost is in fact the product of the prediction probabilities
# MLE did a good job finding parameters that maximize P(D|b0, b1)
table(fit.p > .5)
summary(fit.p[fit.p < .5])
# prediction probabilities of the controls
summary(fit.p[data$Class == 0])
# prediction probabilities of the cases
(fit.p[data$Class == 1])
# all below .5, why is this? 
# model parameters were basically dominated by the data of the controls
# the case data try pulling the parameters to other directions but can only help so much 

# what about adding 100 case
data = rbind(bc0, bc1[sample(1:nrow(bc1), 100), ])
fit <- glm(Class ~ Cell.shape, family = "binomial", data=data)
summary(fit)
(fit.p = predict(fit, type = "response")) %>% hist
table(fit.p > .5)
summary(fit.p[data$Class == 0])
summary(fit.p[data$Class == 1])
# this is much better

# note that model intercept is the over probability w/o any predictors
fit <- glm(Class ~ 1, family = "binomial", data=data)
alpha = coef(fit)[1]
pfun(alpha)
mean(as.numeric(data$Class) - 1)

# note
# logistic regression analysis is essentially finding parameter values that maximize the cost function
# parameter values are further transformed to a 0-1 value (fit.p) by the sigmoid function
# these fit.p values can be interpreted as the event probability based on the given data, although highly imbalanced 
# and fit.p for a given subject will be shifted to the left or right depending on how many total cases in the training data
# this is a probability!
# to use these probabilities for classification is a different thing
# you need to choose a threshold, which is not necessary 0.5
# consequently, different threshold yields different sensitivity and different specificity
# this is a subjective thing

# ROC and AU(ROC) show how specificity and sensitivity go with different thresholds
# and these are metrics to score overall model performance, which itself does not determine the best threshold

# Accuracy is one metric to defind the best threshold, which accounts both sensitivity and specificity

data = rbind(bc0, bc1[sample(1:nrow(bc1), 10), ])
fit <- glm(Class ~ Cell.shape, family = "binomial", data=data)
(fit.p = predict(fit, type = "response")) %>% hist

cutpoints = seq(0, 1, by = 0.01)

accuracy <- sapply(cutpoints, function(cut) {
  fit.y <- ifelse(fit.p > cut, 1, 0)
  fit.y <- factor(fit.y, levels=c(0, 1))
  mean(data$Class == fit.y)
  
})

plot(cutpoints, accuracy, type = "b")

# the curver saturates quickly
# accuracy curve of ATEZO results follow the same pattern, but saturated slower (around 0.5).
