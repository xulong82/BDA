# Decision tree: classification model in the form of a tree structure
# http://www.saedsayad.com/decision_tree.htm

# Learn tree involves selecting input variables and split points on those variables until a suitable tree is constructed

## Algorithm: ID3 (Iterative Dichotomiser 3)

## Entropy of the output variable

entropy <- function(x) { -x * log2(x) - (1-x) * log2(1-x) }
p = seq(from = 0, to = 1, by = 1e-2)
plot(p, entropy(p))

## Entropy of the aoutput variable and a predictor
entropy2 <- function(tbl) tbl

## Information gain: decrease in entropy after a dataset is split on a predictor
## Construcign a decision tree is all about finding predictor that returns the highest information gain
