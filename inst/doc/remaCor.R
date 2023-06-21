## ----setup, echo=FALSE, results="hide"----------------------------------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png",
                      package.startup.message = FALSE,
                      message=FALSE, error=FALSE, warning=TRUE)

## ----rema---------------------------------------------------------------------
library(remaCor)
library(metafor)
library(mvtnorm)
library(clusterGeneration )

# sample size
n = 30

# number of response variables
m = 2

# Error covariance
Sigma = genPositiveDefMat(m)$Sigma

# regression parameters
beta = matrix(0, 1, m)

# covariates
X = matrix(rnorm(n), ncol=1)

# Simulate response variables
Y = X %*% beta + rmvnorm(n, sigma = Sigma)

# Multivariate regression
fit = lm(Y ~ X)

# Correlation between residuals
C = cor(residuals(fit))

# Extract effect sizes and standard errors from model fit
df = lapply(coef(summary(fit)), function(a) 
  data.frame(beta = a["X", 1], se = a["X", 2]))
df = do.call(rbind, df)

# Standard fixed effects meta-analysis
# of independent effects with metafor pacakge
rma( df$beta, sei=df$se, method="FE")

# Standard random effects meta-analysis
# of independent effects with metafor pacakge
rma( df$beta, sei=df$se, method="REML")

# Run fixed effects meta-analysis, assume identity correlation  
# Use Lin-Sullivan method
LS( df$beta, df$se)

# Run fixed effects meta-analysis, accounting for correlation  
# Use Lin-Sullivan method
LS( df$beta, df$se, C)

# Run random effects meta-analysis, assume identity correlation  
RE2C( df$beta, df$se)

# Run random effects meta-analysis, accounting for correlation 
RE2C( df$beta, df$se, C)

## ----out, echo=FALSE----------------------------------------------------------
RE2C( df$beta, df$se, C)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

