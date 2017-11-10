# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1065034/

library(survAUC)
library(survival)

data(pbc)

pbc$event = as.numeric(pbc$status %in% c(1, 2))
pbc$Surv = with(pbc, Surv(time, event))

coxph.age.sex          <- coxph(Surv ~ age + sex, data = pbc)
coxph.age.sex.albumin  <- coxph(Surv ~ age + sex + albumin, data = pbc)

lp.age.sex         <- predict(coxph.age.sex, type = "lp")
lp.age.sex.albumin <- predict(coxph.age.sex.albumin, type = "lp")

TR <- ovarian[1:16,]
TE <- ovarian[17:26,]
train.fit  <- coxph(Surv(futime, fustat) ~ age, x=TRUE, y=TRUE, method="breslow", data=TR)

lpnew <- predict(train.fit, newdata=TE)

GHCI(lpnew)

x = 1:36
beta = seq(from = -1, to = 1, by = 0.1)

beta1 = beta[1]
y = pweibull(x, shape = 1, scale = exp(beta1))
plot(x, y)
