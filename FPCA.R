FPCA = function(simdata,nsamp = 199){
  
  n = length(simdata)
  m = length(simdata[[1]][[5]])
  
  # Extract X and Y data matrices
  X = t(vapply(simdata,function(x) x[[1]],numeric(m)))
  Y = t(vapply(simdata,function(x) x[[2]],numeric(m)))
  
  # Do PCA on X and Y respectively, then project centered X and Y on first n PCs 
  # such that the cumulant percentage is over 99
    prinX     = prcomp(X)
    Mx        = which.min(cumsum(prinX$sdev^2)/sum(prinX$sdev^2)<0.95)
    Xhat      = prinX$x[,1:Mx]

    prinY     = prcomp(Y)
    My        = which.min(cumsum(prinY$sdev^2)/sum(prinY$sdev^2)<0.95)
    Yhat      = prinY$x[,1:My]
    
    set.seed(999)
    dcov.test(Xhat,Yhat, R = nsamp)$p.value
}
