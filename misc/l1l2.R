# E(|x - a|)
# find a to minimize the quantity above

a = seq(-.5, .5, .01)
k = NA

for(i in a) k = c(k, mean(abs(x-i))) # L1
for(i in a) k = c(k, mean((x-i)^2)) # L2

k = k[-1]
plot(a, k)

a[which.min(k)]

mean(x)
median(x)
