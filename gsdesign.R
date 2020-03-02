library(gsDesign)

h0 <- 0.2
h1 <- 0.1
eta <- 0.1
rand.ratio <- 1
tm1 <- 2
tm2 <- 0.5
gamma <- NA
alpha <- 0.05
beta <- 0.1

# 
y1 <- nSurvival(lambda1 = 1/12, lambda2 = 1/24, eta = 0.1, alpha = 0.05, beta = 0.2, Ts = 2, Tr = 1)

#              
y2 <- nSurv(hr = 0.5, hr0 = 1, lambdaC = 1/12, alpha = 0.05, beta = 0.2, eta = 0.1, R = 1, T = 2, minfup = 1)
               
nSurv(lambdaC=log(2)/6, hr=.5, eta=log(2)/40, gamma=6, minfup = 12)
nSurv(lambdaC=log(2)/6, hr=.5, eta=log(2)/40, gamma=3, minfup = 12)

# total events number is derived from the hazard ratio, alpha, beta, and allocation ratio
# from event number to total sample size is a bit complicated

# with censoring rate and hazard rate, you only needs to know the average follow-up time to get the total patient number
# and after that, you can provide either accural rates or accural duration, to derive the other
# but, the average follow-up time is also dependent on the accural duration, or the accural rate
# so, you need to provides one of the two values to get the total patient number

# you can also provide both accural rate and accural duration, hence the sample size, and compute the follow-up time
# but you want to make sure the input sample size is not too big to over-power, and not too small to reach the event number

nSurv(lambdaC=log(2)/6, hr=.5, eta=log(2)/40, gamma= 6, R = 18, minfup = NULL)
