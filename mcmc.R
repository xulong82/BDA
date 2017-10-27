# My take on MCMC
  
# https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
# https://jeremykun.com/2015/04/06/markov-chain-monte-carlo-without-all-the-bullshit/

# Monte Carlo sampling is the cornerstone for the sampling problem, but becomes probihitively inefficient in many applications
# This is where Markov Chain Motne Carlo comes into play, including Metrapolis-Hasting, Gibbs, et al.
# Markov Chain is essentially a random walk on a graph

# MCMC works because of the stationary distribution theorem, aka the fundamental theorem of Markov chains, that states:
# for a very long random walk, the probability that you end at some vertex v is independent of where you started,  
# all of these probabilities taken together is called the stationary distribution of the random walk

# A random walk graph, or Markov chain, corresponds to a transition matrix, from which 
# stationary distribution of the graph vertices is determined and can be computed analytically. 

# In reality, we do not know the transition matrix, while we know the probability of each vertex upon inquiry.
# We certainly do not want to, and often practically impossible, go over each possible value of the variable under study.
# We want to choose a random walk (Markov chain, transition matrix, algorithm) that quickly ouput the stationary distribution
# It is easy to come up with a graph/walk that has the right stationary distribution 
# It is hard to come up with a graph/walk that converges to the stationary distribution fast
# Here comes the Metrapholis-Hasting, which essentially define the transition matrix by using the probability of each vertex upon inquiry

# Metrapolis-Hasting algorithm 
  
# Create data
a <- 0 # intercept
b <- 5 # slope
d <- 5 # standard deviation
n <- 31 # sample number
x <- (-(n-1)/2):((n-1)/2) # independent x-values 
y <-  a + b * x + rnorm(n = n, mean = 0, sd = d) # dependent values according to a + b*x + N(0,sd)
plot(x, y, main = "Test Data")

# Likelihood function
# Model: y = a + b*x + N(0, d) 

likelihood <- function(param) {
  p1 = param[1]
  p2 = param[2]
  p3 = param[3]
  
  prediction = p1 + p2*x
  singlelikelihoods = dnorm(y, mean = prediction, sd = p3, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)
}

## Prior
prior <- function(param){
  p1 = param[1]
  p2 = param[2]
  p3 = param[3]
  p1_prior = dnorm(p1, mean = 0, sd = 5, log = T)
  p2_prior = dunif(p2, min = 0, max = 10, log = T)
  p3_prior = dunif(p3, min = 0, max = 10, log = T)
  return(p1_prior + p2_prior + p3_prior)
}

## Posterior 
posterior <- function(param){
  return (likelihood(param) + prior(param))
}

# Metropolis algorithm
# 1. Starting at a random parameter value
# 2. Choosing a new parameter value close to the old value based on some probability density that is called the proposal function
# 3. Jumping to this new point with a probability p(new)/p(old), where p is the target function, and p>1 means jumping as well

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1, 3))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = rnorm(3, mean = chain[i,], sd = c(0.5,0.1,0.3))
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    
    if (runif(1) < probab) {
      chain[i+1,] = proposal
    } else {
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

startvalue = c(0, 5, 10)
chain = run_metropolis_MCMC(startvalue, 1e4)

burnIn = 5e3
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

# Plotting

par(mfrow = c(3,1))
hist(chain[-(1:burnIn), 1], main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnIn), 1]))
abline(v = a, col="red" )
hist(chain[-(1:burnIn), 2], main="Posterior of b", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn), 2]))
abline(v = b, col="red" )
hist(chain[-(1:burnIn), 3], main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn), 3]) )
abline(v = d, col="red" )

# for comparison:
summary(lm(y ~ x))
