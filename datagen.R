## Loading packages
library(mvtnorm)
library(wavethresh)

#Covariance Structure #Suppose K=8, k=1,...,K
  # Rho is the correlation structure of xi's and zeta's
  # betaX and betaY are defined as shrinking speed of eigenvalues
  # n is the size of data

###############################################################################
#                                   Normal                                    #
###############################################################################

# Eigenvalues decay polynomially
coefGenPoly = function(Rho, betaX, betaY, n) {
  K = nrow(Rho)
  Sigma11 = diag(4 * seq(1, K) ^ (-betaX))
  Sigma22 = diag(4 * seq(1, K) ^ (-betaY))
  Sigma12 = sqrt(Sigma11) %*% Rho %*% sqrt(Sigma22)
  Sigma = rbind(cbind(Sigma11, Sigma12),
                cbind(t(Sigma12), Sigma22))
  Coef = rmvnorm(n, sigma = Sigma) # First K columns are xi for X, second K columns are zeta for Y
  varX = sum(diag(Sigma11))
  varY = sum(diag(Sigma22))
  list(Coef,varX,varY)
}

# Eigenvalues increase polynomially
coefGenPolyRvrs = function(Rho, betaX, betaY, n) {
  K = nrow(Rho)
  Sigma11 = diag(4 * seq(K, 1, -1) ^ (-betaX))
  Sigma22 = diag(4 * seq(K, 1, -1) ^ (-betaY))
  Sigma12 = sqrt(Sigma11) %*% Rho %*% sqrt(Sigma22)
  Sigma = rbind(cbind(Sigma11, Sigma12),
                cbind(t(Sigma12), Sigma22))
  Coef = rmvnorm(n, sigma = Sigma) # First K columns are xi for X, second K columns are zeta for Y
  varX = sum(diag(Sigma11))
  varY = sum(diag(Sigma22))
  list(Coef,varX,varY)
}

###############################################################################
#                      Uncorrelated but Dependent                             #
###############################################################################

coefGenUncorPoly = function(betaX, n, K) {
  sigma = diag(4 * seq(1, K) ^ (-betaX))
  coefX = rmvnorm(n, sigma = sigma)
  coefY = vapply(1:K,
                 function(k) {
                   coefX[, k] ^ 2 - sigma[k, k]
                 }, numeric(n))
  Coef = cbind(coefX, coefY) # First K columns are xi for X, second K columns are zeta for Y
  varX = sum(diag(sigma))
  varY = 2 * sum(diag(sigma)^2) # chi-square variance
  list(Coef,varX,varY)
}

coefGenUncorMix = function(betaX, betaY, n, K, cut) {
  sigma = diag(c(4 * seq(1, K) ^ (-betaX), 4 * seq(1, cut)^(-betaY)))
  coefXY1 = rmvnorm(n, sigma = sigma)
  coefY2  = vapply((cut + 1):K,
                   function(k) {
                     coefXY1[, k] ^ 2 - sigma[k, k]
                   }, numeric(n))
  Coef = cbind(coefXY1, coefY2) # First K columns are xi for X, second K columns are zeta for Y
  varX = sum(diag(sigma)[1:K])
  varY = sum(diag(sigma)[(K+1):(K+cut)]) + 2 * sum((diag(sigma)[(cut+1):K])^2)
  list(Coef,varX,varY)
}

###############################################################################
#                     Simulation Data Generating function                     #
###############################################################################
simDataGen = function(basisX, basisY, Coef, Time, SNR) {
  #############################################################################
  # Input:
  ##    basisX and basisY are basis functions of X's and Y's respectively
  ##    Coef[[1]] is a matrix of coefficients for basis function of X's and Y's
  ##    Coef[[2]] is sum of eigenvalues for X
  ##    Coef[[3]] is sum of eigenvalues for Y
  ##    Time is a list of sampling time for subject i=1,...,n
  ##    SNR is signal to noise ratio defined by Karhunen-Loeve expansion

  # Output:
  ## A list of n sublist which contains X_ij, Y_ij and time t_ij,
  
  #############################################################################
  
  n = length(Time) #actually number of subject
  K = ncol(Coef[[1]]) / 2
  xi = Coef[[1]][, 1:K]
  coefsX = lapply(1:n, function(x) t(xi[x, ]))
  zeta = Coef[[1]][, (K + 1):(2 * K)]
  coefsY = lapply(1:n, function(x) t(zeta[x, ]))
                  
  mapply(
    FUN = function(X,Y,t){
        # X_noi = X + rnorm(length(X),sd = sqrt(0.1 * var(X)))
        # Y_noi = Y + rnorm(length(Y),sd = sqrt(0.1 * var(Y)))
        noise_X = rnorm(length(X),sd = sqrt(Coef[[2]]/SNR))
        noise_Y = rnorm(length(Y),sd = sqrt(Coef[[3]]/SNR))
        X_noi = X + noise_X
        Y_noi = Y + noise_Y
        X_w = threshold(wd(X_noi,family = "DaubExPhase"),policy = "BayesThresh")
        Y_w = threshold(wd(Y_noi,family = "DaubExPhase"),policy = "BayesThresh")
        list(
             wr(X_w),
             wr(Y_w),
             X_w$D,
             Y_w$D,
             t,
             wd(X_noi,family = "DaubExPhase")$D,
             wd(Y_noi,family = "DaubExPhase")$D
             #wd(noise_X,family = "DaubExPhase")$D,
             #wd(noise_Y,family = "DaubExPhase")$D
            )
    },
    mapply(
      FUN = function(b, c) {
        as.vector(c %*% b)
      },
      lapply(Time, basisX),
      coefsX,
      SIMPLIFY = FALSE
    ),
    mapply(
      FUN = function(b, c) {
        as.vector(c %*% b)
      },
      lapply(Time, basisY),
      coefsY,
      SIMPLIFY = FALSE
    ),
    Time,
    SIMPLIFY = FALSE
  )
}
