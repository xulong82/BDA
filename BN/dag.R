library(abn)
library(rstan)

# discover the best fitting DAG for Iris dataset
irisData <- data.frame(model.matrix(~Sepal.Length + Sepal.Width + Petal.Length + Petal.Width + Species - 1, data=iris))
irisData$Speciessetosa = NULL

# binomial nodes to factors
irisData$Speciesversicolor <- factor(irisData$Speciesversicolor)
irisData$Speciesvirginica <- factor(irisData$Speciesvirginica)

# ban arcs between generated dummy variables
banMat <- matrix(0, nrow = ncol(irisData), ncol = ncol(irisData))
banMat[5, 6] = banMat[6, 5] = 1
rownames(banMat) <- colnames(banMat) <- names(irisData)

# distributions for each node
dists <- as.list(c(rep("gaussian", 4), rep("binomial", 2)))
names(dists) = names(irisData)

# build auxiliary cache
cache <- buildscorecache(data.df = irisData, data.dists = dists, dag.banned = banMat, max.parents = 5)

# find the best fitting DAG
DAG <- mostprobable(cache)

# fit the model, compute the marginals and make a nice plottable graph
abnModel <- fitabn(dag.m = DAG, data.df = irisData, data.dists = dists, create.graph = T, compute.fixed = T)

# plot the DAG
plot(abnModel$graph)

# check the posterior modes
abnModel$modes

# Stan!!!

# standardize the data.
irisData[1:4] <- scale(irisData[1:4])

# also, we want use it as a classifier for the species, so let's build a model with training and test data
train <- sample(150, 100, replace=F)
trainSet <- irisData[train,]
testSet <- irisData[-train,]

# fit and sample from the model
dataStan <- list(Ntest = dim(testSet)[1], Ntrain=dim(trainSet)[1],
                 sLength1 = trainSet$Sepal.Length, sWidth1 = trainSet$Sepal.Width,
                 pLength1 = trainSet$Petal.Length, pWidth1 = trainSet$Petal.Width,
                 sVersi = as.numeric(trainSet$Speciesversicolor)-1, sVirgi = as.numeric(trainSet$Speciesvirginica)-1,
                 sLength = testSet$Sepal.Length, sWidth = testSet$Sepal.Width,
                 pLength = testSet$Petal.Length, pWidth = testSet$Petal.Width)
modelStan <- stan_model(file='dag.stan')
fitStan <- sampling(modelStan, data = dataStan, iter=500, warmup=250)

# extract pars. We want to compare the most probable class to the true class.
pars <- extract(fitStan)

probVirgi <- apply((pars$virgi/(pars$virgi+pars$notVirgi)), 2, mean)
probVersi <- apply((pars$versi/(pars$versi+pars$notVersi)), 2, mean)
probSetosa <- 1- (probVirgi+probVersi)

postProbs <- data.frame(probSetosa, probVersi, probVirgi)

table(apply(postProbs, 1, function(x) x[1]>x[2] & x[1]>x[3]), 
      (testSet$Speciesvirginica == 0) & (testSet$Speciesversicolor == 0))
table(apply(postProbs, 1, function(x) x[2]>x[1] & x[2]>x[3]), 
      testSet$Speciesversicolor)
table(apply(postProbs, 1, function(x) x[3]>x[1] & x[3]>x[2]), 
      testSet$Speciesvirginica)
