library(truncnorm)
library(mvtnorm)

Compute.residuals = function(theta.hat,data.uni.sim)
{
shape = exp(theta.hat[1])
scale = exp(theta.hat[2])

obs.time = data.uni.sim$obs.time
delta = data.uni.sim$delta

N=dim(data.uni)[1]

indices1 = (delta==1) ### observed
indices2 = (delta==0) ### right-censored
indices3 = (delta==-1) ### left-censored

residuals = rep(-999,N)
residuals[indices1] = qnorm(pweibull(obs.time[indices1],shape=shape,scale=scale,lower.tail=FALSE))
residuals[indices2] = etruncnorm(b=qnorm(pweibull(obs.time[indices2],shape=shape,scale=scale,lower.tail=FALSE)))
residuals[indices3] = etruncnorm(a=qnorm(pweibull(obs.time[indices3],shape=shape,scale=scale,lower.tail=FALSE)))

return(residuals)
}


Compute.test.stat.vect = function(theta.hat,data.uni.sim,indices.biv,data.biv,h1=1/2)
{
residuals = Compute.residuals(theta.hat,data.uni.sim)
resid1 = residuals[indices.biv[,1]]
resid2 = residuals[indices.biv[,2]]
k1 = data.biv$kinship.ff
k2 = data.biv$kinship.fm+data.biv$kinship.mf+data.biv$kinship.mm
k = h1*k1 + (1-h1)*k2
result = resid1*resid2*k
result.fam = aggr(result,data.biv$fam.id)
result.fam
}

Compute.test.stat = function(theta.hat,data.uni.sim,indices.biv,data.biv)
{
R1 = genD(Compute.test.stat.vect,theta.hat,data.uni.sim=data.uni.sim,indices.biv=indices.biv,data.biv=data.biv)

U = R1$f0
C = apply((R1$D)[,1:2],2,mean)

R2 = genD(Log.Lik.Vect.H0,theta.hat,data.uni.sim=data.uni.sim)$D

first.deriv = R2[,1:2]
second.deriv = R2[,3:5]

B = solve(matrix(apply(second.deriv,2,mean)[c(1,2,2,3)],ncol=2))

corr.term = (t(C)%*%B%*%t(first.deriv))[1,]

stat = sum(U)
var.stat = sum((U-corr.term)^2)
Z.obs = stat/sqrt(var.stat)

return(Z.obs)
}



##############
## Minimum p  ##
##############


prepare.test.pmin = function(theta.hat,data.uni.sim,indices.biv,data.biv)
{
R1 = genD(Compute.test.stat.vect,theta.hat,data.uni.sim=data.uni.sim,indices.biv=indices.biv,data.biv=data.biv,h1=1)
R2 = genD(Compute.test.stat.vect,theta.hat,data.uni.sim=data.uni.sim,indices.biv=indices.biv,data.biv=data.biv,h1=0)

U1 = R1$f0
C1 = apply((R1$D)[,1:2],2,mean)

U2 = R2$f0
C2 = apply((R2$D)[,1:2],2,mean)

R3 = genD(Log.Lik.Vect.H0,theta.hat,data.uni.sim=data.uni.sim)$D

first.deriv = R3[,1:2]
second.deriv = R3[,3:5]

B = solve(matrix(apply(second.deriv,2,mean)[c(1,2,2,3)],ncol=2))

W1 = (t(C1)%*%B%*%t(first.deriv))[1,]
W2 = (t(C2)%*%B%*%t(first.deriv))[1,]

return(data.frame(U1=U1,U2=U2,W1=W1,W2=W2))
}

test.pmin=function(theta.hat,data.uni.sim,indices.biv,data.biv,seq.h=seq(0,1,0.1))
{
UW = prepare.test.pmin(theta.hat,data.uni.sim,indices.biv,data.biv)
U1=UW$U1
U2=UW$U2
W1=UW$W1
W2=UW$W2

U = c(sum(U1),sum(U2))
Var.U = matrix(c(sum((U1-W1)^2),sum((U1-W1)*(U2-W2)),sum((U1-W1)*(U2-W2)),sum((U2-W2)^2)),ncol=2)

matr.V = cbind(seq.h,1-seq.h)
V = matr.V%*%U
Var.V = matr.V%*%Var.U%*%t(matr.V)

matr.Z = diag(1/sqrt(diag(Var.V)))
Z = matr.Z%*%V
Var.Z = matr.Z%*%Var.V%*%t(matr.Z)

pmin=min(1-pnorm(Z))
pvalue.min = 1-pmvnorm(upper=rep(max(Z),length(seq.h)),sigma=Var.Z)

return(c(pmin,pvalue.min,1-pnorm(Z)))
}

