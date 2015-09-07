library(rstan)
library(Biobase)
library(genefilter)
library(dplyr)
library(tissuesGeneExpression)
library(limma)
library(rafalib)

rm(list = ls())
setwd("~/Dropbox/Stan")
source("../X/function.R")

load("gene_nonp.rdt")
y <- as.matrix(log2(gene_nonp + 1))
g <- factor(gsub("(M|W).*", "\\1", colnames(gene_nonp)), levels = c("W", "M"))
sds <- rowSds(y)

# estimate prior distribution
ebayes <- rstan::stan_model(file = "./scaled_inverse_chi_square_y.stan", model_name="ebayes")
dat <- list(n = length(sds), y = sds) 
(fit <- sampling(ebayes, dat, warmup = 3e2, iter = 2e3, chains = 3))
sample <- as.data.frame(fit)[c("nu", "tau")]
(mode <- apply(sample, 2, function(z) { dens <- density(z); dens$x[which.max(dens$y)] }))

# test prior estimation
ebayes <- rstan::stan_model(file = "./scaled_inverse_chi_square.stan", model_name="ebayes")
dat <- list(nu = 14, tau = 0.45) # estimation
dat <- list(nu = 30, tau = 0.45) # adjustment: better than estimation, why?
(fit <- sampling(ebayes, dat, warmup = 3e2, iter = length(sds), chains = 3))
sample <- as.data.frame(fit)$sigma

plot(density(sds), xlim = c(0, 1))
lines(density(sample), xlim = c(0, 1), col = "red")

# plug in prior, estimate posterior
ebayes <- rstan::stan_model(file = "./ebayes.stan", model_name="ebayes")
dat <- list(N = 6, g = as.numeric(g), nu = 30, tau = 0.45)

dat$y <- y[1, ]
(fit <- sampling(ebayes, dat, warmup = 3e2, iter = 6e2, chains = 3))
(lm0 <- lm(y[1, ] ~ g) %>% summary)

sample <- as.data.frame(fit)[c("alpha", "beta", "sigma")]
(mode <- apply(sample, 2, function(z) { dens <- density(z); dens$x[which.max(dens$y)] }))
(ci95 <- summary(fit)$summary[c("alpha", "beta", "sigma"), c("2.5%", "97.5%")])

ci95 <- apply(y, 1, function(x) {
  dat1 = within(dat, y <- x)
  fit <- sampling(ebayes, dat1, warmup = 3e2, iter = 6e2, chains = 3)
  summary(fit)$summary["beta", c("2.5%", "97.5%")]
}) %>% t

gene_select <- rownames(ci95)[rowMin(ci95) > 0.0]
gene_select <- rownames(ci95)[rowMin(ci95) > 0.2]
gene_select <- rownames(ci95)[rowMin(ci95) > 0.3]
gene_select <- rownames(ci95)[rowMin(ci95) > 0.4]
gene_select <- rownames(ci95)[rowMin(ci95) > 0.5]
gene_select <- rownames(ci95)[rowMin(ci95) > 1.0]

gk <- mmGK(gene_select)
data.frame(KEGG = gk$KEGG$Term[1:20], BP = gk$GO$BP$Term[1:20], MF = gk$GO$BP$Term[1:20])

# --- limma method (moderated t-statistics) ---
data("tissuesGeneExpression")
# MLE by optim() does not work well!
log.scaled.f <- function(x, df, par) -sum(log(pf(x / par[1], df, par[2])))
optim(par = c(0.1, 5), log.scaled.f, x = biosds^2, df = 11)
x = optim(par = c(0.1, 5), log.scaled.f, x = biosds^2, df = 11)$par
tmp=hist(biosds^2, xlab="sd",ylab="density",freq=FALSE,nc=800,xlim=c(0,0.4))
sds = seq(0,0.4,len=100)
dd = df(sds/x[1], 11, x[2])
k=sum(tmp$density)/sum(dd) # normalizing constant 
lines(sds,dd*k,type="l",col=3,lwd=2)

estimates=fitFDist(biosds^2,11)
theoretical<- sqrt(qf((seq(0,999)+0.5)/1000,11,estimates$df2)*estimates$scale)
observed <- biosds
qqplot(theoretical,observed)

###### 
# model: hierarchical linear model
# prior on variance: scaled F distribution
# posterior function? 

y = e[, which(tissue %in% c("endometrium", "cerebellum"))]
g = as.factor(tissue[tissue %in% c("endometrium", "cerebellum")])
var = rowVars(y)

lm.fit <- lmFit(y, design=model.matrix(~ g))
ebayes.fit = eBayes(lm.fit)

df.prior = fdist.fit$df2
df.total = df + df.prior -1
df.total = df.residual + out$df.prior

var.prior = fdist.fit$scale
var.post = (df*var + df.prior*var.prior) / df.total #	squeeze posterior variances

out$t = coefficients / stdev.unscaled / sqrt(out$s2.post)
out$p.value = 2*pt(-abs(out$t),df=df.total)

coefficients <- lm.fit$coefficients
stdev.unscaled <- lm.fit$stdev.unscaled # what is this?
df.residual <- lm.fit$df.residual
apply(y[1:10, ], 1, function(x) lm(x ~ g)$coef)
coefficients[1:10, ]

sampleSD = ebayes.fit$sigma # residual standard error
posteriorSD = sqrt(ebayes.fit$s2.post) # sigma as prior to s2.post

plot(sampleSD, posteriorSD)
abline(0, 1)

ebayes <- function(fit, proportion=0.01, stdev.coef.lim=c(0.1,4)) {
  coefficients <- fit$coefficients
  stdev.unscaled <- fit$stdev.unscaled
  sigma <- fit$sigma
  df.residual <- fit$df.residual
  
  out <- squeezeVar(sigma^2, df.residual)
  out$s2.prior <- out$var.prior
  out$s2.post <- out$var.post
  out$var.prior <- out$var.post <- NULL
  df.total <- df.residual + out$df.prior
  out$t <- coefficients / stdev.unscaled / sqrt(out$s2.post)
  out$p.value <- 2*pt(-abs(out$t),df=df.total)
}

squeezeVar <- function(var, df)
{
  out <- fitFDist(var, df1 = df)
  out$var.prior <- out$scale
  out$df.prior <- out$df2
  out$df2 <- out$scale <- NULL
  df.total <- df + out$df.prior
  out$var.post <- (df*var + out$df.prior*out$var.prior) / df.total
  out
}

fitFDist <- function(x, df1) {
  z <- log(x) #	Better to work on with log(F)
  e <- z-digamma(df1/2)+log(df1/2)
  emean <- mean(e)
  evar <- mean(n/(n-1)*(e-emean)^2-trigamma(df1/2))
  if(evar > 0) {
    df2 <- 2*trigammaInverse(evar)
    s20 <- exp(emean+digamma(df2/2)-log(df2/2))
  } else {
    df2 <- Inf
    s20 <- exp(emean)
  }
  list(scale=s20,df2=df2)
}

N = 30
X.ns = mean(bwt.nonsmoke)
sd.ns = sd(bwt.nonsmoke)
X.s = mean(bwt.smoke)
sd.s = sd(bwt.smoke)
sd.diff = sqrt(sd.ns^2/N+sd.s^2/N)
tval = (X.ns - X.s)/sd.diff
