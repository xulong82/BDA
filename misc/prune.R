library(igraph)
library(MASS)

# Objective: identify minimal cutsets that disconnect all edges in a graph
# We approach the problem based on topologies of graphs

# Algorithm
# 1. Find a node (N0) with maximal connectivity
# 2. Remove N0 and update the graph
# 3. Repeat 1-2 until no edges are left (complete fragmentation)

# --- playground ---

toy = graph_from_literal( A--B, C--D, E--F, G--H, I, J, K )
toy = make_graph("Bull")
plot(toy)

# set vertex name if needed
V(toy)
toy = set_vertex_attr(toy, "name", value = letters[1:length(V(toy))])
plot(toy)

# find cliques
cliques(toy, min = 3)

# total edges
ecount(toy)

# adjacency number of each vertex
g.adj <- get.adjacency(toy) # matrix
v.top <- which.max(rowSums(as.matrix(g.adj)))
names(v.top)

# delete the top vertex
toy = delete.vertices(toy, names(v.top)) 
plot(toy)

# --- the full algorithm ---

toy = make_graph("Bull")
plot(toy)

g.del = NULL

while(ecount(toy) > 0) {
  g.adj <- get.adjacency(toy)
  v.top <- which.max(rowSums(as.matrix(g.adj)))
  g.del <- c(g.del, names(v.top))
  toy = delete.vertices(toy, names(v.top)) 
}

plot(toy)

# I think this algorithm works. 
# Please check if it breaks on certain graph structures.

# Now, let's review Leon's VIF idea
# VIF (variance inflation factor) quantifies collinearity of the predictors in regression models  

# Toy data (100 by 5), with collinearity

vif.dt = as.data.frame(matrix(rnorm(5e2, 0, 1), nrow = 1e2))
vif.dt$V2 = vif.dt$V1 + rnorm(1e2, 0, .3)
vif.dt$V3 = vif.dt$V1 + rnorm(1e2, 0, .6)
vif.dt$V4 = vif.dt$V1 + rnorm(1e2, 0, .9)
cor(vif.dt)

# compute the VIF of each variables

(vif <- sapply(names(vif.dt), function(x) {
  toy.dt = vif.dt[c(x, setdiff(names(vif.dt), x))]
  names(toy.dt) = paste0("V", 1:5)
  toy.lm = summary(lm(V1 ~ V2 + V3 + V4 + V5, data = toy.dt))
  1 / ( 1 - toy.lm$r.squared )
}))

# a convenient way
(vif = diag(MASS::ginv(cor(vif.dt))))
names(vif) = names(vif.dt)

# keep the variable with top VIF
# as expected the top variable is V1, giving the simulation procedure
(vif.top = which.max(vif))

# we want to remove variables that are highly correlated with the top variable
# to do so, we will re-compute VIF after taking off the top variable, and compute the VIF drop

vif.new = diag(MASS::ginv(cor(vif.dt[-vif.top])))
names(vif.new) = names(vif.dt)[-vif.top]

# the new VIF will drop because the regression models drop an important predictor
(vif.drop = vif[-vif.top] - vif.new)

# as expected, VIF of V1 drops the most, but is this significant?
# we can obtain empirical null distribution of the VIF drops, by repeatedly permuting the top variable 
# this way, we repeatedly plug in a different top variable, while preserving the relationship of the remaining variables

vif.drop.null <- t(replicate(1e3, {
  vif.dt[, vif.top] = sample(vif.dt[, vif.top])
  vif = diag(MASS::ginv(cor(vif.dt)))
  vif.new = diag(MASS::ginv(cor(vif.dt[-vif.top])))
  vif[-vif.top] - vif.new
}))

# set significant cut at the .95 percentile
(vif.drop.cut <- apply(vif.drop.null, 2, function(x) quantile(x, .95)))
vif.drop - vif.drop.cut

# V2, V3, V4 will be excluded
# Final outcomes will be V1 and V5

# Leon actually permute all the remaining variables
vif.drop.null.leon <- t(replicate(1e3, {
  toshuffle = setdiff(1:ncol(vif.dt), vif.top)
  vif.dt[, toshuffle] = sapply(toshuffle, function(x) sample(vif.dt[, x])) 
  vif.new = diag(MASS::ginv(cor(vif.dt[-vif.top])))
  vif[-vif.top] - vif.new
}))

# and set a single significant cut for all remaining variables
(vif.drop.cut.leon <- quantile(c(vif.drop.null.leon), .95))
vif.drop - vif.drop.cut.leon

# none will be excluded 

# Remarks: 
# This VIF idea looks unnecessarily cumbersome. 
# Plus, the null hypothesis needs to be clarified in order to decide the exact permutation procedure
