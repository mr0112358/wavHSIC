dCov_c = function(simdata,nsamp = 199){
  
  n = length(simdata)
  m = length(simdata[[1]][[5]])
  
  # Extract X and Y data matrices
  X = t(vapply(simdata,function(x) x[[1]],numeric(m)))
  Y = t(vapply(simdata,function(x) x[[2]],numeric(m)))
  
  set.seed(999)
  dcorT.test(X,Y)$p.value
}
