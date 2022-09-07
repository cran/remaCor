## ----setup, echo=FALSE, results="hide"----------------------------------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png",
                      package.startup.message = FALSE,
                      message=FALSE, error=FALSE, warning=TRUE)

## ----rema---------------------------------------------------------------------
library(remaCor)
library(metafor)

# Generate effects
library(mvtnorm)
library(clusterGeneration )

n = 4
Sigma = cov2cor(genPositiveDefMat(n)$Sigma)
beta = t(rmvnorm(1, rep(0, n), Sigma))
stders = rep(.1, n)

# Standard fixed effects meta-analysis
# of independent effects with metafor pacakge
rma( beta, sei=stders, method="FE")

# Standard random effects meta-analysis
# of independent effects with metafor pacakge
rma( beta, sei=stders, method="REML")

# Run fixed effects meta-analysis, assume identity correlation  
# Use Lin-Sullivan method
LS( beta, stders)

# Run fixed effects meta-analysis, accounting for correlation  
# Use Lin-Sullivan method
LS( beta, stders, Sigma)

# Run random effects meta-analysis, assume identity correlation  
RE2C( beta, stders)

# Run random effects meta-analysis, accounting for correlation 
RE2C( beta, stders, Sigma)

## ----out, echo=FALSE----------------------------------------------------------
RE2C( beta, stders, Sigma)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

