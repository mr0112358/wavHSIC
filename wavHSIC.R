wavHSIC = function(simdata, ratio = 1, nsamp = 199, alphaX = 2.8, alphaY = 2.8){
  # p=q in this version
  
  n = length(simdata)
  m = length(simdata[[1]][[5]])
  
  # For wavelets
  nlevel = log2(m)-1
  
  # Extract data
  ## Denoised coeffs
  X_dn    = t(vapply(simdata, function(x) {rev(x[[3]])}, numeric(m-1)))
  Y_dn    = t(vapply(simdata, function(x) {rev(x[[4]])}, numeric(m-1)))
  X_n     = t(vapply(simdata, function(x) {rev(x[[6]])}, numeric(m-1)))
  Y_n     = t(vapply(simdata, function(x) {rev(x[[7]])}, numeric(m-1)))

# X_dn, Y_dn are coefficients matrix of X and Y, each row is a sample curve, scales are ordered from 1 to J
# Return log2 of dVar(X) and dVar(Y) in different scales
  log2Dvar = lapply(0:nlevel, function(k){
          distX_dn = dist(x = as.matrix(X_dn[,seq(2^k,2^(k+1)-1)]), method = "minkowski", p = 2, upper = T)
          distY_dn = dist(x = as.matrix(Y_dn[,seq(2^k,2^(k+1)-1)]), method = "minkowski", p = 2, upper = T)
          distX_n = dist(x = as.matrix(X_n[,seq(2^k,2^(k+1)-1)] - X_dn[,seq(2^k,2^(k+1)-1)]), method = "minkowski", p = 2, upper = T)
          distY_n = dist(x = as.matrix(Y_n[,seq(2^k,2^(k+1)-1)] - Y_dn[,seq(2^k,2^(k+1)-1)]), method = "minkowski", p = 2, upper = T)
          result_dn = dcovU_stats(as.matrix(distX_dn),as.matrix(distY_dn))
          result_n = dcovU_stats(as.matrix(distX_n),as.matrix(distY_n))
          return(list(distX_dn,distY_dn,c(log2(result_dn[3]), log2(result_dn[4])), c(as.integer(result_n[3] < ratio * result_dn[3]), as.integer(result_n[4] < ratio * result_dn[4]))))
      })
    
  cutJ = rowSums(vapply(1:(nlevel+1), function(k) return(log2Dvar[[k]][[4]]), numeric(2)))
  J_X = cutJ[1]
  J_Y = cutJ[2]
    
  #a_X = (log2Dvar[[2]][[3]][1] - log2Dvar[[J_X]][[3]][1])/(J_X-2)/2
  #a_Y = (log2Dvar[[2]][[3]][2] - log2Dvar[[J_Y]][[3]][2])/(J_Y-2)/2
  
  gammaX = sapply(1:J_X, function(x) log2Dvar[[x]][[3]][1])
  gammaY = sapply(1:J_Y, function(x) log2Dvar[[x]][[3]][2])
  JX = 1:J_X
  JY = 1:J_Y
  a_X = -coef(lm(gammaX ~ JX))[2]/2
  a_Y = -coef(lm(gammaY ~ JY))[2]/2
  
  # Calculate Distance Matrices
  distX = sqrt(Reduce('+', lapply(1:(nlevel+1), function(k){
      (log2Dvar[[k]][[1]]*2^(k*a_X))^2
  })))
  distY = sqrt(Reduce('+', lapply(1:(nlevel+1), function(k){
      (log2Dvar[[k]][[2]]*2^(k*a_Y))^2
  })))
  set.seed(999)
  c(pvalue = dcov.test(distX,distY,index=1, R = nsamp)$p.value, beta_X = a_X, beta_Y = a_Y)
}

