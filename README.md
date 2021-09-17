# Replication Code for Entropic Latent Variable Integration via Simulation [(Schennach, 2014)](https://doi.org/10.3982/ECTA9748)
## Summary
This repository is my attempt at reproducing the estimation procedure from Susanne Schennach's paper titled Entropic Latent Variable Integration via Simulation (ELVIS) published in Econometrica (2014). The original code is written in GAUSS and translated to the R programming language. The estimator is a simulation-based algorithm which uses an MCMC procedure to integrate out unobserved heterogeneity from models defined by unconditional moment restrictions using a posterior distribution found by entropy maximization. The files contained here reproduce estimates from an interval-valued data regression, whose parameters are known to be partially identified. 

There are a number of computational improvements that need to be made before the estimator can be generalized to other models. The main difficulty in implementation is knowing a priori whether the model considered is point or set (partially) identified. While in theory, the estimator proposed in this paper adapts to either of these cases, estimators for models that are point identified are computationally much simpler than those that are set-identified, for example using [Chernozhukov, Hong, and Tamer (2007)](https://doi.org/10.1111/j.1468-0262.2007.00794.x). Future updates to this repository will consider both of these cases.




