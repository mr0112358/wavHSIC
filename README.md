# Code Example for "A Wavelet-Based Independence Test for Functional Data with an Application to MEG Functional Connectivity"

## Method functions
- [wavHSIC.R](./wavHSIC.R) wavHSIC, the proposed method
- [FPCA.R](./FPCA.R) dCov on FPC scores with 95\% variation
- [dcov_c.R](./dcov_c.R) Distance correlation t-test of independence in high dimension ([Szekley and Rizzo, 2013](https://www.sciencedirect.com/science/article/pii/S0047259X13000262))
- [Pearson.R](./Pearson.R) Pearson's test ([He et al. 2012, sec 2.4.3, steps 1 & 2](https://www.sciencedirect.com/science/article/abs/pii/S0197458011005744))
- [gtemp.R](./gtemp.R) Global Temporal Correlation ([Zhou et al., 2018](https://www.sciencedirect.com/science/article/pii/S0047259X1730341X))
- [dnm.R](./dnm.R) Dynamical Correlation ([Dubin and Mueller, 2005](https://amstat.tandfonline.com/doi/abs/10.1198/016214504000001989))
- [KMSZ.R](./KMSZ.R) Testing for lack of dependence in the functional linear model ([Kokoszka et al., 2008](https://onlinelibrary.wiley.com/doi/abs/10.1002/cjs.5550360203))

## Data Generation
Noised Data given SNR: [datagen.R](./datagen.R)

## Simulation for different settings
- An example is given in [example.R](./example.R).
