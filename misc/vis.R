rm(list = ls())
file_path <- "http://www.sthda.com/sthda/RDoc/data/housetasks.txt"
housetasks <- read.delim(file_path, row.names = 1)

library(gplots)
library(graphics)
library(corrplot)

# Graphical display of contengency tables

dt <- as.table(as.matrix(housetasks))
gplots::balloonplot(t(dt), main ="housetasks", xlab ="", ylab="", label = FALSE)
corrplot::corrplot(as.matrix(housetasks), is.cor = FALSE)
graphics::mosaicplot(dt, shade = TRUE, las=2, main = "housetasks")

# Chi-squared test

(chisq <- chisq.test(housetasks))
chisq$observed
chisq$expected
chisq$residuals

# Pearson residuals, aka standardized residuals, r $\frac{o-e}{e}$, reflects how much the observed values deviate from the expected.
corrplot::corrplot(chisq$residuals, is.cor = FALSE)

contrib <- 100 * chisq$residuals^2/chisq$statistic
round(contrib, 3)

corrplot::corrplot(contrib, is.cor = FALSE)

# it can be seen that the most contributing cells to the X2 are Wife/Laundry (7.74%), Husband/Repairs (21.9%), Jointly/Holidays (12.44%);
# these cells contribute about 47.06% to the total X2 score and thus account for most of the difference between expected and observed values;
# visual interpretation may be complex when the contingency table is very large; 
# in this case, the contribution of one cell to the total X2 score becomes a useful way of establishing the nature of dependency.
