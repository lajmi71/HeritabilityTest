library(numDeriv)
source("Functions.1.R")

Log.Lik.Vect.H0 = function(theta,data.uni.sim)
{
shape = exp(theta[1])
scale = exp(theta[2])

N=nrow(data.uni.sim)

LogLik = rep(-999,N)

indices1 = (data.uni.sim$ind.id == 1)  ## proband
indices2 = (data.uni.sim$ind.id != 1) & (data.uni.sim$delta == 1) ## non proband with observed event
indices3 = (data.uni.sim$ind.id != 1) & (data.uni.sim$delta == 0) ## non proband with right-censored event
indices4 = (data.uni.sim$ind.id != 1) & (data.uni.sim$delta == -1) ## non proband with left-censored event

LogLik[indices1] = dweibull(data.uni.sim$obs.time[indices1],shape=shape,scale=scale,log=TRUE)-pweibull(data.uni.sim$assessmentAge[indices1],shape=shape,scale=scale,log.p=TRUE)
LogLik[indices2] = dweibull(data.uni.sim$obs.time[indices2],shape=shape,scale=scale,log=TRUE)
LogLik[indices3] = pweibull(data.uni.sim$obs.time[indices3],shape=shape,scale=scale,lower.tail=FALSE,log.p=TRUE)
LogLik[indices4] = pweibull(data.uni.sim$obs.time[indices4],shape=shape,scale=scale,log.p=TRUE)

aggr(A=LogLik,B=data.uni.sim$fam.id)
}

Neg.Log.Lik.H0 = function(theta,data.uni.sim) {-sum(Log.Lik.Vect.H0(theta,data.uni.sim))}

Estim.H0 = function(theta,data.uni.sim) {optim(theta,Neg.Log.Lik.H0,data.uni.sim=data.uni.sim)$par}

Estim.var.theta.hat.H0 = function(theta.hat,data.uni.sim)
{
R=genD(Log.Lik.Vect.H0,theta.hat,data.uni.sim=data.uni.sim)$D

N = dim(R)[1]

first.deriv = R[,1:2]
second.deriv = R[,3:5]

A=(t(first.deriv)%*%first.deriv)/N
B = solve(matrix(apply(second.deriv,2,mean)[c(1,2,2,3)],ncol=2))

B%*%A%*%t(B)/N
}