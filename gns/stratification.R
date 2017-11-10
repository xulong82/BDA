library(dplyr)
library(data.table)
library(glmnet)
library(REFSfs)
library(BMS)
library(mclust)
library(nnet)

# ? cluster the subjects

# fs results on 1124 patients :: 132 networks :: 100 samples :: 2 arms
load("~/roche34/otf2/fs1pfs.rdt") # starmaster
load("~/Projects/roche34/forward/fs1pfs.rdt")
nrow(myfs1) == 1124 * 132 * 100 * 2

# summarize the fs samples per patient per condition
fs = data.table(myfs1)
fs = fs[, c("condition", "fixedDataRow", "network", "output", "Variable_Output")]
fs$output = log2(fs$output + 1)
fs.stats = fs[, list(output.avg = mean(output), output.sd = sd(output)), by = "condition,fixedDataRow"]

fs.stats.r <- fs.stats[fs.stats$condition == "R", ] 
fs.stats.g <- fs.stats[fs.stats$condition == "G", ] 
all(fs.stats.r$fixedDataRow == fs.stats.g$fixedDataRow)

# differential effects per patient per network
fs.diff = data.frame(sample = fs.stats.r$fixedDataRow, diff = fs.stats.g$output.avg - fs.stats.r$output.avg)
fs.diff$pval = pval

# --------------------------------------
# ? any difference by taking only 100 samples to derive the effects per patient per arm
fs.diff2 = sapply(c("R", "G"), function(x1) {
  fs.x1 = fs[fs$condition == x1, ]
  sapply(1:1124, function(x2) { cat(x1, x2, "\n")
    fs.x2 = fs.x1[fs.x1$fixedDataRow == x2, ]
    fs.x3 = fs.x2[sample(1:nrow(fs.x2), size = 100, replace = F), "output"]
    mean(fs.x3$output)
  })
})

fs.diff = data.frame(sample = 1:1124, diff = log2(fs.diff2[, "G"] + 1) - log2(fs.diff2[, "R"] + 1))
fs.diff$pval = pval
# --------------------------------------

with(fs.diff, plot(diff, pval, xlab = "Differential Effects (G versus R)"))
with(fs.diff, plot(abs(diff), pval))
with(fs.diff, hist(diff, n = 20))
with(fs.diff, hist(pval, n = 20))
with(fs.diff, hist(-log10(pval), n = 20))

# gaussian mixture decomposition
mclust = Mclust(fs.diff$diff)
mclust = Mclust(as.matrix(fs.diff))
plot(mclust, what = "classification")
plot(mclust, what = "density") # overlay with hist(fs.diff$diff)
plot(mclust, what = "BIC") # to understand it

load("~/roche34/otf2/alldata.rdt") # training data
ens = fsReadModel("/users/xwang/roche34/otf2/m1.new/NetFragFile.txt.protobuf") # network ensemble
edges = fsEdgeFrequencies(ens, freqThreshold = 0) # edges
edges.ooi = filter(edges, output == "pfs_AVAL")
edges.ins = as.character(edges.ooi$input)

mylrt <- sapply(edges.ins, function(x) {
  lrt.df = data.frame(de = fs.diff$diff, biomarker = refsdf.c1[, x])
  loglk0 = logLik(lm(de ~ 1, data = lrt.df))
  loglk1 = logLik(lm(de ~ 1 + biomarker, data = lrt.df))
  loglk1 - loglk0
})

mylrt = data.frame(biomarker = edges.ins, LRT = mylrt)
mylrt$P.value = 1 - pchisq(2 * mylrt$LRT, df = 1)

mylrt = mylrt[order(mylrt$P.value), ]
mylrt$P.adj = p.adjust(mylrt$P.value, method = "bonferroni", n = ncol(refsdf.c1))  
rownames(mylrt) = NULL

save(mylrt, file = "~/roche34/fs/mylrt.rdt")
# scp starmaster:~/roche34/fs/mylrt.rdt ~/Projects/roche34/forward

load("~/Projects/roche34/forward/mylrt.rdt")
datatable(mylrt)

edges.ins = as.character(mylrt$biomarker)
mylrt2 <- sapply(edges.ins, function(x) {
  lrt.df = data.frame(de = mclust$classification, biomarker = refsdf.c1[, x])
  loglk0 = logLik(multinom(de ~ 1, data = lrt.df))
  loglk1 = logLik(multinom(de ~ 1 + biomarker, data = lrt.df))
  loglk1 - loglk0
})

mylrt2 = data.frame(biomarker = edges.ins, LRT = mylrt2)
mylrt2$P.value = 1 - pchisq(2 * mylrt2$LRT, df = 1)

mylrt2 = mylrt2[order(mylrt2$P.value), ]
mylrt2$P.adj = p.adjust(mylrt2$P.value, method = "bonferroni", n = ncol(refsdf.c1))  
rownames(mylrt2) = NULL

datatable(mylrt2)

plot(mylrt$LRT, mylrt2$LRT)
temp = merge(mylrt, mylrt2, by = "biomarker")

plot(-log10(temp$P.value.x), -log10(temp$P.value.y), xlab = "-log10(P) Continuous", ylab = "-log10(P) Discrete")
abline(0, 1)

# test multiple biomarkers by regularized regression
mylasso <- glmnet(y=y.1, x= x.1, family="gaussian")

# multinomial regression using cluster identity
multinom(y ~ x1 + x2, df1)

# hierarchical clustering
mydist = dist(fs.out$diff, method = "euclidean")
myhclust = hclust(mydist)
plot(myhclust, labels = F)
# subpopulation identification by different cut values
mycutree = cutree(myhclust, k = 2:50) # cut trees to sub-groups
# make design matrix by using cutree results
mydesign = apply(mycutree, 2, function(x) {
  ids = unique(x)
  mat = matrix(0, nrow = length(x), ncol = length(ids))
  for ( id in ids) mat[x == id, id] = 1 
  mat[, -1] # redundancy
})
mybic = sapply(mydesign, function(x) BIC(lm(fs.out$diff ~ x)))
plot(names(mybic), mybic)
# bic formula
loglik = logLik(lm(fs.out$diff ~ mydesign[[1]]))
log(nrow(fs.out)) * 3 - 2 * loglik

# biomarker identification
# regularized (lasso, ipredict) multinomial regression:
# cluster_identity ~ candidate_predictors

# LRT method per network 
fs = data.table(myfs1)
fs = fs[, c("condition", "fixedDataRow", "network", "output", "Variable_Output")]
fs = fs[, list(output.avg = mean(output), output.med = median(output)), by = "condition,fixedDataRow,network"]

fs.g1 <- fs[fs$condition == "R", ] 
fs.g2 <- fs[fs$condition == "G", ] 

all(fs.g1$network == fs.g2$network)
all(fs.g1$fixedDataRow == fs.g2$fixedDataRow)
  
# differential effects per patient per network by LRT
fs.out = mutate(fs.g1, diff = log2(fs.g1$output.avg + 1) - log2(fs.g2$output.avg + 1)) 
fs.out$condition = NULL
  
fs.out.net = fs.out[fs.out$network == 1, ]
  
edges.ins = as.character(edges.ooi$input)
mylrt <- sapply(edges.ins, function(x) {
  lrt.df = data.frame(de = fs.out.net$diff, biomarker = refsdf.c1[, x])
  loglk0 = logLik(lm(de ~ 1, data = lrt.df))
  loglk1 = logLik(lm(de ~ 1 + biomarker, data = lrt.df))
  loglk1 - loglk0
})

mylrt = data.frame(biomarker = edges.ins, lrt = mylrt)
mylrt$P.value = 1 - pchisq(2 * mylrt$lrt, df = 1)

# differential effects per patient per network by PIP
# initialize results data.frame
load("~/roche34/otf2/alldata.rdt") # training data
ens = fsReadModel("/users/xwang/roche34/otf2/m1.new/NetFragFile.txt.protobuf") # network ensemble
bayesian.res = data.frame(Predictor = numeric(0), PostIncluProb = numeric(0), PostExpectEffect = numeric(0), Network = numeric(0))
  
ioi = "asl.treat.char_ARMCD"
ooi = "pfs_AVAL"

alldata = refsdf.c1

# retrieve causal vars for each network
lapply(1:128, function(i) { 
  net1 = fsSubsetEnsemble(ens, i)
  frag = fsGetFrags(net1, "pfs_AVAL")
  frag$input %>% unlist
})

for (i in 1:128) { cat("network:", i, "\n")
  net1 = fsSubsetEnsemble(ens, i)
# vars = fsCausalVars(net1, ooi, cutoff = 0, maxpath = -1)
  vars = fsCausalVars(net1, ooi, cutoff = 0, maxpath = 1)
  vars = vars[! vars %in% ioi]
    
  if (length(vars) == 0) {
    next
  } else if (length(vars) == 1) {
    data0 = alldata[vars]
    data0$diff = fs.out[fs.out$network==i,]$diff
      
    glm0 = lm(data0$diff ~ 1)
    glm1 = lm(data0$diff ~ 1 + data0[, vars])
    bic0 = BIC(glm0)/2.0
    bic1 = BIC(glm1)/2.0
      
    pip = exp(-bic1)/(exp(-bic1)+exp(-bic0))
    post.mean = pip *coef(glm1)[2][[1]]

    bayesian.res = rbind(bayesian.res, c(vars, pip, post.mean, i))
  } else {
    data0 = alldata[vars]
    data0$diff = fs.out[fs.out$network==i,]$diff
      
    bms1 = bms(diff ~ ., data = data0, mprior = "uniform", g = "UIP", user.int = F, mcmc="enumerate")
    bms1.co = data.frame(coef(bms1))
      
    temp = data.frame(cbind(rownames(bms1.co), bms1.co$PIP, bms1.co$Post.Mean, i))
    bayesian.res = rbind(bayesian.res, temp)
  }
}

names(bayesian.res) = c("Predictor", "PostIncluProb", "PostExpectEffect", "Network")
  
# Produce E[PosteriorInclusionProb|D] for each variable by summing the PosteriorInclusionProb for each variable, 
# and dividing by number of networks
  
posterior_network_results$PosteriorInclusionProb  = as.numeric(as.character(posterior_network_results$PosteriorInclusionProb))
posterior_network_results$PosteriorExpectedEffect = as.numeric(as.character(posterior_network_results$PosteriorExpectedEffect))
  
posterior_results = aggregate(cbind(PosteriorInclusionProb, PosteriorExpectedEffect) ~ causalPredictor, FUN = sum, data=posterior_network_results, na.action = NULL)
  
posterior_results$PosteriorExpectedEffect = posterior_results$PosteriorExpectedEffect/as.numeric(number_networks)
posterior_results$PosteriorInclusionProb  = posterior_results$PosteriorInclusionProb/as.numeric(number_networks)
posterior_results                         = posterior_results[order(-posterior_results$PosteriorInclusionProb),]
  
posterior_results
