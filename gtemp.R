gtemp = function(simdata,nsamp = 199){
    
n = length(simdata)
m = length(simdata[[1]][[5]])

# Extract Data
X = t(vapply(simdata,function(x) x[[1]],numeric(m)))
Y = t(vapply(simdata,function(x) x[[2]],numeric(m)))

# Here we use |rhot| rather than rhot because we just want to test dependency
rhot = vapply(1:m, function(j) abs(cor(X[,j],Y[,j])), numeric(1))
rho = mean(rhot)

# Permutation test
gmap = split(replicate(nsamp,sample.int(n)),rep(1:nsamp, each = n))
allgtempcor = vapply(gmap, function(perm){
  mean(vapply(1:m, function(j) abs(cor(X[,j],Y[perm,j])), numeric(1)))
},numeric(1))
pvalue = mean(allgtempcor > rho)
}
