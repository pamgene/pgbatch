library(pgBatch)
library(R6)
X = matrix(nrow = 100, ncol = 20, data = rnorm(20*100))
bx = as.factor(c( rep(1,10) ,rep(2,10) ))
grp =
cMod = pgCombat$new()
cMod = cMod$fit(X, bx,ref.batch = "1")
Xc = cMod$apply(X, bx)
Xc = cMod$Xc
Xcc = sva::ComBat(X, bx, ref.batch = "1")

print(all(round(Xc,8) - round(Xcc,8) == 0))
