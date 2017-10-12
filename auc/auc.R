library(pROC)

# auc: auc of a classifier is equivalent to the probability that the classifier 
# will rank a randomly chosen positive instance higher than a randomly chosen negative instance

# Titanic data

df = read.csv("~/github/statistics/auc/titanic3.csv")

df$pclass <- as.factor(df$pclass)
df$survived <- as.factor(df$survived)
df$sex <- as.factor(df$sex)

df <- df[,c("survived","pclass","sex","age","sibsp","parch")]
df <- df[complete.cases(df),]

N <- nrow(df)
size=10 # leave-10-out CV
  
df <- df[sample(N), ] # shuffle
  
num <- floor(N/size)
rest <- N - num * size
ncv <- cumsum(c(rep(size,num), rest))
  
predictions <- data.frame(survived = df$survived, pred = NA)
  
for(n in ncv) {
  v <- rep(TRUE, N)
  v[(n-size+1):n] <- FALSE
  
  lr <- glm(survived ~ ., data = df[v,], family = binomial(logit))
  predictions[!v,"pred"] <- predict(lr, newdata=df[!v,], type="response") 
  # type = "resposne" gives predicted probabilities
  # default is the linear predictors (log-odds of the probabilities)
}

v <- rep(NA, nrow(predictions))
v <- ifelse(predictions$pred >= threshold & predictions$survived == 1, "TP", v)
v <- ifelse(predictions$pred >= threshold & predictions$survived == 0, "FP", v)
v <- ifelse(predictions$pred < threshold & predictions$survived == 1, "FN", v)
v <- ifelse(predictions$pred < threshold & predictions$survived == 0, "TN", v)

predictions$pred_type <- v

# distribution of predicted survival probabilities as violin and jitter points
ggplot(data= predictions, aes(x=survived, y=pred)) + 
  geom_violin(fill=rgb(1,1,1,alpha=0.6), color=NA) + 
  geom_jitter(aes(color=pred_type), alpha=0.6) +
  scale_color_discrete(name = "type") +
  labs(title=sprintf("Threshold at %.2f", threshold))

pROC::auc(predictions$survived, predictions$pred)
