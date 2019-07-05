aggr = function(A,B)
{
res=NULL
for(b in unique(B))
{
res=c(res,sum(A[B==b]))
}
return(res)
}