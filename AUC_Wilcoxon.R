# AUC meets u-statistics

# AUC is directly connected to the Mann-Whitney u-statistic
# AUC = U / N1 / N2

category <- c(1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0)
prediction <- rev(seq_along(category))
prediction[9:10] <- mean(prediction[9:10])

library('pROC')
(official_auc <- auc(roc(category, prediction)))

idx <- as.logical(category)
g1 <- prediction[idx]
g2 <- prediction[!idx]
u = wilcox.test(g1, g2)$statistic

u/(length(g1) * length(g2))
