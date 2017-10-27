library(dplyr)
library(data.table)
library(glmnet)
library(REFSfs)
library(BMS)

load("~/roche34/otf2/alldata.rdt")

ens = fsReadModel("/users/xwang/roche34/otf2/m1.new/NetFragFile.txt.protobuf") # network ensemble
n.nets = fsNumNets(ens)

edges = fsEdgeFrequencies(ens, freqThreshold = 0) 
terms = fsTermFrequencies(ens, incParamStats=T, freqThreshold = 0)

ooi = "pfs_AVAL"
edges.ooi = filter(edges, output == ooi)

load("~/roche34/otf2/fs1pfs.rdt")
load("~/Projects/roche34/forward/fs1pfs.rdt")
# forward simulation results on 1124 patients :: 132 networks :: 100 samples :: 2 arms
nrow(myfs1) == 1124 * 132 * 100 * 2

# step 1: take differential effects results by forward simulation
fs = data.table(myfs1)
fs = fs[, c("condition", "fixedDataRow", "network", "output", "Variable_Output")]
fs = fs[, list(output.avg = mean(output), output.med = median(output)), by = "condition,fixedDataRow"]

fs.g1 <- fs[fs$condition == "R", ] 
fs.g2 <- fs[fs$condition == "G", ] 

all(fs.g1$fixedDataRow == fs.g2$fixedDataRow)

# differential effects per patient per network
fs.out = mutate(fs.g1, diff = log2(fs.g1$output.avg + 1) - log2(fs.g2$output.avg + 1)) 
# fs.out = mutate(fs.g1, diff = log2(fs.g1$output.med + 1) - log2(fs.g2$output.med + 1)) 
fs.out$condition = NULL

edges.ins = as.character(edges.ooi$input)
mylrt <- sapply(edges.ins, function(x) {
  lrt.df = data.frame(de = fs.out$diff, biomarker = refsdf.c1[, x])
  loglk0 = logLik(lm(de ~ 1, data = lrt.df))
  loglk1 = logLik(lm(de ~ 1 + biomarker, data = lrt.df))
  loglk1 - loglk0
})

mylrt = data.frame(biomarker = edges.ins, LRT = mylrt)
mylrt$P.value = 1 - pchisq(2 * mylrt$LRT, df = 1)

mylrt = mylrt[order(mylrt$P.value), ]
mylrt$P.adj = p.adjust(mylrt$P.value, method = "bonferroni", n = ncol(refsdf.c1))  

# regularized regression
mylasso <- glmnet(y=y.1, x= x.1, family="gaussian")
  
# ------------
# step 2: hierarchical clustering
mydist = dist(fs.out$diff, method = "euclidean")
myhclust = hclust(mydist)
plot(myhclust, labels = F)

# step 3: identify optimal cluster number
# we will get different number of clusters by different cut values
# effect ~ a_variable_of_cluster_identity
# essentially the number of binary dummy variables 
# assess fittness by BIC penalization and use results to identify optimal cluster numbers

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
#-------------

# Revised steps 2-3: Gaussian mixture model decomposition
# Line 3900 in refsfsutils package

library(mclust)

x1 = rnorm(100, 1, 1)
x2 = rnorm(100, 3, 1)

x3 = rnorm(100, 1, 1)
x4 = rnorm(100, 3, 1)

y1 = c(x1, x2)
y2 = c(x3, x4)

mclust = Mclust(cbind(y1, y2))
mclustbic = mclustBIC(cbind(y1, y2))

plot(mclust, what = "classification")
plot(mclustbic)
 
?Mclust
?mclustBIC

# step 4: biomarker identification
# regularized (lasso, ipredict) multinomial regression:
# cluster_identity ~ candidate_predictors
# cvglmnet function in the code

# LRT method per network 
fs = data.table(myfs1)
fs = fs[, c("condition", "fixedDataRow", "network", "output", "Variable_Output")]
fs = fs[, list(output.avg = mean(output), output.med = median(output)), by = "condition,fixedDataRow,network"]

fs.g1 <- fs[fs$condition == "R", ] 
fs.g2 <- fs[fs$condition == "G", ] 

all(fs.g1$network == fs.g2$network)
all(fs.g1$fixedDataRow == fs.g2$fixedDataRow)
  
# differential effects per patient per network
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

# dan method
# initialize results data.frame
bayesian.res = data.frame(Predictor = numeric(0), PostIncluProb = numeric(0), PostExpectEffect = numeric(0), Network = numeric(0))
  
ioi = "asl.treat.char_ARMCD"
ooi = "pfs_AVAL"

alldata = refsdf.c1

# retrieve causal vars for each network
lapply(1:n.nets, function(i) { 
  net1 = fsSubsetEnsemble(ens, i)
  frag = fsGetFrags(net1, "pfs_AVAL")
  frag$input %>% unlist
})

for (i in 1:n.nets) { cat("network:", i, "\n")
  net1 = fsSubsetEnsemble(ens, i)
  vars = fsCausalVars(net1, ooi, cutoff = 0, maxpath = -1)
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
      
    bms1 = bms(diff ~ ., data = data1, mprior = "uniform", g = "UIP", user.int = F, mcmc="enumerate")
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
