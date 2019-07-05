
compute.biv.indices=function(data.uni,kin1,kin2)
{
indice.ind1=NULL
indice.ind2=NULL
k1 = NULL
k2 = NULL
compteur=0
for(i in unique(data.uni$fam.id))
{
compteur = compteur+1
kin1.fam = kin1[[compteur]]
kin2.fam = kin2[[compteur]]
ind.id.fam = data.uni$ind.id[data.uni$fam.id==i]
len = length(ind.id.fam)
for (j in 1:(len-1))
{
for(k in (j+1):len)
{
indice.ind1 = c(indice.ind1,which((data.uni$fam.id==i)&(data.uni$ind.id==ind.id.fam[j])))
indice.ind2 = c(indice.ind2,which((data.uni$fam.id==i)&(data.uni$ind.id==ind.id.fam[k])))
k1 = c(k1,kin1.fam[j,k])
k2 = c(k2,kin2.fam[j,k])
}
}
}
cbind(indice.ind1,indice.ind2,k1,k2)
}
