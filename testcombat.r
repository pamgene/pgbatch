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
cMod = cMod$fit(X, bx, mean.only = FALSE)
Xc = cMod$apply(X, bx)
Xc1 = cMod$Xc
Xcc = sva::ComBat(X, bx, mean.only  = FALSE)

print(all(round(Xc,8) - round(Xcc,8) == 0))

dfs = df %>% group_by(RunID, ID) %>% dplyr::summarize(mval =mean(value))
Xm = acast(dfs, ID ~ RunID, value.var = "mval")
bxm = acast(dfs, ID~ RunID, value.var = "RunID")[1,]
mMod = pgCombat$new()
mMod = cMod$fit(Xm, bxm, mean.only = TRUE)
Xcm = cMod$apply(Xm, bxm)
Xccm =  sva::ComBat(Xm, bxm, mean.only  = TRUE)
"C:/Users/rdwijn/Documents/210-000/210-226 PamStationDX/Control/Ref samples/Verification/RefveR"
