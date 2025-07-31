
prob = c(0.75, 0.95)

varper=c(50,90)
kper <- varper/100
varfun="ns"
ppred<-sort(c(1:99,97.5))/100
argvarm<-list(fun=varfun, knots=kper)
argvar <- list(x=ppred,fun=varfun, knots=kper, Bound=c(0.01,0.99))
bvar <- do.call(onebasis,argvar)

nlag <- 4; breaksl<-1
xlag <- 0:nlag
arglagm<-list(fun="strata", breaks=breaksl)
cen<-"m"