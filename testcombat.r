library(pgBatch)
library(R6)
library(reshape2)
library(dplyr)

df = combatapp::combat_testdf

X = acast(df, rowSeq~colSeq, value.var = "value")
bx = acast(df, rowSeq~colSeq, value.var = "RunID")[1,]
bx = factor(bx)

#X = matrix(nrow = 100, ncol = 22, data = rnorm(22*100))
#dimnames(X) = list(rowSeq = 1:dim(X)[1], colSeq = 1:dim(X)[2])
#bx = as.factor(c( rep(1,7) ,rep(2,7), rep(3,8) ))


cMod = pgCombat$new()
cMod = cMod$fit(X, bx, mean.only = TRUE,ref.batch = "1")
Xc = cMod$apply(X, bx)
Xc1 = cMod$Xc
Xcc = sva::ComBat(X, bx, mean.only  = TRUE, ref.batch = "1")

print(all(round(Xc,8) - round(Xcc,8) == 0))



