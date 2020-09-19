dnm = function(simdata,nsamp = 199){
  
  n = length(simdata)
  m = length(simdata[[1]][[5]])
  
  # Extract Data
  X = t(vapply(simdata,function(x) x[[1]],numeric(m)))
  Y = t(vapply(simdata,function(x) x[[2]],numeric(m)))
  
  # Coefficients in Simpson's formula
  SimpVec = rep(1/m,times = m)
  
  # Standardize X's and Y's
  Xcup = X - matrix(rep(colMeans(X),each = n),ncol = m)
  Mx = Xcup %*% SimpVec
  diffXcup = Xcup - matrix(rep(Mx,times=m),ncol = m)
  Xstar = matrix(rep(1/sqrt(diffXcup^2 %*% SimpVec),times=m),ncol = m) * diffXcup
  
  Ycup = Y - matrix(rep(colMeans(Y),each = n),ncol = m)
  My = Ycup %*% SimpVec
  diffYcup = Ycup - matrix(rep(My,times=m),ncol = m)
  Ystar = matrix(rep(1/sqrt(diffYcup^2 %*% SimpVec),times=m),ncol = m) * diffYcup
  
  # Here we use |rho| rather than rho because we just want to test dependency
  dnmcor = mean((Xstar * Ystar) %*% SimpVec)
  
  # Permutation test
  gmap = split(replicate(nsamp,sample.int(n)),rep(1:nsamp, each = n))
  alldnmcor = vapply(gmap, function(perm){
    mean((Xstar * Ystar[perm,]) %*% SimpVec)
  },numeric(1))
  pvalue = mean(abs(alldnmcor)>=abs(dnmcor))
}
