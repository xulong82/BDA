library(survival)
library(survAUC)
library(Hmisc)

# y1: predicted response
# y2: observed response

y1 = rnorm(100, 10, 1)
y2 = rnorm(100, 10, 1)

concordance <- function(y1, y2)
{
  s <- 0 # score
  n <- 0 # number of pairs
  for(i in seq(along=y1))
  {
    for(j in seq(along=y2))
    {
      if(i != j)
      {
        if(y1[i] > y1[j])
        {
          s <- s + (y2[i] > y2[j]) + 0.5 * (y2[i] == y2[j])
          n <- n + 1
        }
      }
    }
  }
  s / n
}

concordance(y1, y2)
Hmisc::rcorr.cens(y1, y2) # match with Hmisc method

cens = sample(c(0, 0, 0, 1), size = 100, replace = T)
mysurv = Surv(y1, cens)

Hmisc::rcorr.cens(y2, mysurv)
survAUC::UnoC(mysurv, mysurv, y2)
