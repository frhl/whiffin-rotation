
# a model where orfs and prt are correlated, but
# only lens causally affect protein levels
k <- 1000
lens <- rnorm(k, 1)
orfs <- rnorm(k, lens*0.5)
prt <- rnorm(k, lens*0.5)


summary(lm(prt ~ orfs))
summary(lm(prt ~ orfs + lens))



#u5 len ->  prt
# orf -> prt



