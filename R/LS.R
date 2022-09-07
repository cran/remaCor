# Adapted from RE2C from http://software.buhmhan.com/RE2C/index.php
# by Gabriel Hoffman
# December 19, 2019

#' Fixed effect meta-analysis for correlated test statistics
#'
#' Fixed effect meta-analysis for correlated test statistics using the Lin-Sullivan method.
#'
#' @param beta regression coefficients from each analysis
#' @param stders standard errors corresponding to betas
#' @param cor correlation matrix between of test statistics.  Default considers uncorrelated test statistics 
#'
#' @details Perform fixed effect meta-analysis for correlated test statistics using method of Lin and Sullivan (2009).  By default, correlation is set to identity matrix to for independent test statistics.  
#'
#' This method requires the correlation matrix to be symmatric positive definite (SPD).  If this condition is not satisfied, results will be NA.  If the matrix is not SPD, there is likely an issue with how it was generated. 
#'
#' However, evaluating the correlation between observations that are not pairwise complete can give correlation matricies that are not SPD.  In this case, consider running \code{Matrix::nearPD( x, corr=TRUE)} to produce the nearest SPD matrix to the input. 
#'
#' @references{
#'   \insertRef{lin2009meta}{remaCor}
#' }
#'
#' @return
#' Return values:
#' \itemize{
#' \item{\code{beta}: }{effect size}
#' \item{\code{se}: }{effect size standard error}
#' \item{\code{p}: }{p-value}
#'}
#'
#' @examples
#' # Generate effects
#' library(mvtnorm)
#' library(clusterGeneration )
#' 
#' n = 4
#' Sigma = cov2cor(genPositiveDefMat(n)$Sigma)
#' beta = t(rmvnorm(1, rep(0, n), Sigma))
#' stders = rep(.1, n)
#' 
#' # Run fixed effects meta-analysis, 
#' # assume identity correlation  
#' LS( beta, stders)
#' 
#' # Run random effects meta-analysis,
#' # assume identity correlation  
#' RE2C( beta, stders)
#' 
#' # Run fixed effects meta-analysis, 
#' # account for correlation 
#' LS( beta, stders, Sigma)
#' 
#' # Run random effects meta-analysis,
#' # account for correlation 
#' RE2C( beta, stders, Sigma)
#'
#' @import stats
#' @export
LS <- function(beta, stders, cor=diag(1, length(beta)) ){

   # check arguments
   if( length(beta) != length(stders) ){
      stop("Number of test statistics and standard errors must be the same")
   }
   if( nrow(cor) != ncol(cor) ){      
      stop("Correlation matrix must be square")
   }
   if( length(beta) != nrow(cor) ){      
      stop("Number of test statistics and rows of correlation matrix must be the same")
   }
   if( ! is.numeric(beta) ){      
      stop("beta must be numeric")
   }
   if( ! is.numeric(stders) ){      
      stop("stders must be numeric")
   }
   if( any(stders <= 0) ){      
      stop("All values in stders must be positive")
   }

   ## conventional FE approach  
   V <- diag(stders) %*% cor %*% diag(stders)

   # invert V, if not valid then warn and return empty results
   # Vinv <- solve(V)
   Vinv <- tryCatch( solve(V),
    error = function(e){
      warning("At least 1 eigen-value of correlation matrix (cor) is negative,\nso the matrix is not a valid (i.e. positive definite) correlation matrix.\nConsider using Matrix::nearPD().")
      NA
    }
  )  
  if( length(Vinv) == 1 ){
   if( is.na(Vinv) )
      return( data.frame(beta = NA, se = NA, p = NA) )
  }

   ones <- matrix(rep(1,length(beta)),nrow=1)
   
   newx <- (ones %*% Vinv %*% beta) / (ones %*% Vinv %*% t(ones))
   newv <- 1 / (ones %*% Vinv %*% t(ones))
   
   if( newv > 0 ){    
      newstd <- sqrt(newv)  
      newp <- pchisq(newx*newx/newv, 1, lower.tail=FALSE)
   }else{
      warning("At least 1 eigen-value of correlation matrix (cor) is negative,\nso the matrix is not a valid (i.e. positive definite) correlation matrix.\nConsider using Matrix::nearPD().")
      newstd = NA
      newp = NA
   }

   data.frame(beta = newx, se = newstd, p = newp)
}




