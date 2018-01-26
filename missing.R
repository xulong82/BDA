x <- rnorm(100)
y <- 2*x + 3

x2 <- rnorm(30)
y2 <- rep(0, 30)

yy <- c(y, y2)
xx <- c(x, x2)

i <- c(rep(1, 100), rep(0, 30))

zz = i * xx

z1 = lm(yy ~ 0 + i + i:xx)
z2 = lm(yy ~ 0 + i + zz)
z3 = lm(yy[1:100] ~ 0 + i[1:100] + xx[1:100])

# How to treat missing values?
# Cons and Pros

# 1. delete the samples
# 2. imputation
# 3. categorize the variables
