# Load packages
pack<-installed.packages()
packages<-pack[,1] 
if(!is.element("energy", packages)){install.packages("energy")}
if(!is.element("wavethresh", packages)){install.packages("wavethresh")}
if(!is.element("mvtnorm", packages)){install.packages("mvtnorm")}
library(energy)

# Load Parameters
n = 50 # Sample size
m = 256 # Number of time points
K= 16 # Number of basis 
SNR = 8 # Signal to noise ratio
albe = c(1.05,1.2) # var(eta) = k^{-1.05}, var(zeta) = k^{-1.2}
setID = 2 # can be 1,2 or 3

Setting = list(Setting1 = list("Indep", "Normal", diag(rep(0,16))),
               Setting2 = list("PartDepdt", "Normal", diag(c(rep(0,8),rep(0.6,8)))),
               Setting3 = list("PartUncor", "Uncor", 8))

# Load source files
source("datagen.R") #data generation
source("Pearson.R")
source("dnm.R")
source("gtemp.R")
source("dcov_c.R")
source("FPCA.R")
source("KMSZ.R")
source("wavHSIC.R")

# Set time grid
Time = rep(list(0:(m-1) / m), times = n)

# Basis functions
sqrt2 = sqrt(2)
X_basis = function(t) {
  sqrt2 * t(matrix(unlist(lapply(seq(2, K, 2), function(k) {
    c(cos(k * pi * t), sin(k * pi * t))
  })),
  ncol = K))
}
Y_basis = function(t) {
  sqrt2 * t(matrix(unlist(lapply(seq(2, K, 2), function(k) {
    c(cos(k * pi * (t + 0.2)), sin(k * pi * (t + 0.2)))
  })),
  ncol = K))
}

# Generate simulation data for different cases
if(setID==3){
  cut = Setting[[setID]][[3]]
  betaX = albe[1]
  betaY = albe[2]
  Coef = coefGenUncorMix(betaX, betaY, n, K, cut)
  set.seed(11) # can be changed
  simdata = simDataGen(X_basis, Y_basis, Coef, Time, SNR)
}else{
  Rho = Setting[[setID]][[3]]
  betaX = albe[1]
  betaY = albe[2]
  Coef = coefGenPoly(Rho, betaX, betaY, n)
  set.seed(11) # can be changed
  simdata = simDataGen(X_basis, Y_basis, Coef, Time, SNR)
}

# Output
Pearson(simdata) #return p-value
dnm(simdata)     #return p-value
gtemp(simdata)   #return p-value
dCov_c(simdata)  #return p-value
FPCA(simdata)    #return p-value
KMSZ(simdata)    #return p-value
wavHSIC(simdata) #return p-value, beta_X and beta_Y

