source("Estim.H0.1.R")
source("Sim.1.R")
source("Create.kinship.1.R")
source("Tests.1.R")
source("indices.1.R")

fam.id = data.uni$fam.id  ### family indentifiant
ind.id = data.uni$ind.id  ### individual indentifiant
assessmentAge = data.uni$assessmentAge ### Age at assessment for each individual

unique.fam.id = unique(fam.id)

#####################################
## Parameters to generate the data ##
#####################################

### Shape and scale of marginal Weibull distribution 

shape = 3.5   
scale = 140
theta=c(log(shape),log(scale))

### Heritability parameters

h1=0.7
h2=0


data.uni.sim = Generate.data(shape,scale,kin1,kin2,h1,h2,fam.id,ind.id,assessmentAge,unique.fam.id)  ##### Generate data set
indices.biv = compute.biv.indices(data.uni.sim,kin1,kin2)

theta.hat = Estim.H0(theta,data.uni.sim) ##### Estimate the parameters under the null hypothesis 

pvalue1 = 1-pnorm(Compute.test.stat(theta.hat,data.uni.sim,indices.biv,data.biv))  ##### Testing heritability under model 1
pvalue2 = test.pmin(theta.hat,data.uni.sim,indices.biv,data.biv,seq.h=seq(0,1,0.1))[2] ##### Testing heritability under model 2