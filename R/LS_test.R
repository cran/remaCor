

# LS.emp.resid = function(beta, stders, P, nu, n.mc.samples=1e4, seed=1){

# 	lambda = corpcor::estimate.lambda(P, verbose=FALSE)

# 	S = (1-lambda)*cor(P) + lambda*diag(diag(cor(P)), ncol(P))

# 	# compute Lin-Sullivan test-statistic
# 	res = LS(beta, stders, S)

# 	# so that crossprod(A) is the covariance of beta 
# 	A = (scale(P) / sqrt(n)) %*% diag(stders)

# 	# Compute empirical p-value using a gamma fit to 
# 	# Monte Carlo samples from the null distribution
# 	res$p = .LS_empirical_pvalue(with(res, (beta/se)^2), A, nu,stders, lambda, n.mc.samples, seed )

# 	res
# }


# #' @importFrom EnvStats egamma
# .LS_empirical_pvalue = function(LS.stat, A, nu, stders, lambda, n.mc.samples=1e4, seed=1){

# 	stopifnot( all(nu > 1) )

# 	# sample stats assuming finite sample size
# 	stats = get_stat_samples(A, nu, stders, lambda, n.mc.samples, seed)

# 	# estimate gamma parameters
# 	fit.g = egamma(stats)

# 	# compute p-value
# 	pgamma(LS.stat, 
# 		shape = fit.g$parameters[1], 
# 		scale=fit.g$parameters[2], lower.tail=FALSE)
# }

# get_stat_samples = function(A, nu, stders, lambda, n.mc.samples, seed){

# 	if( exists(".Random.seed") ){
# 		old <- .Random.seed
# 		on.exit({.Random.seed <<- old})
# 	}
# 	set.seed(seed)

# 	ch = chol(crossprod(A))

# 	sapply(seq(n.mc.samples), function(i) get_sample_new(ch, nu, stders, lambda, seed))
# }

# get_sample_new = function(ch, nu, stders, lambda, seed){

# 	# dcmp = svd(A)
# 	# 1) Sample beta from MVN
# 	# beta_i = rmvnorm( 1, rep(0, ncol(A)), crossprod(A))
# 	# beta_i = matrix(rnorm(ncol(ch)), ncol=ncol(ch)) %*% ch
# 	beta_i = Rfast::matrnorm(1, ncol(ch)) %*% ch
# 	beta_i = c(beta_i)

# 	# 2) Sample variance from Wishart
# 	# this is the cholesky of the covariance matrix
# 	# C_i_chol = remaCor:::rwishart_chol(nu, ch) / sqrt(nu)
# 	# V_i = crossprod(C_i_chol)
	
# 	# P_prime = rmvnorm( n, rep(0, ncol(ch)), crossprod(ch)) / sqrt(nu)
# 	# P_prime = matrix(rnorm(n*ncol(ch)), ncol=ncol(ch)) %*% ch / sqrt(nu)
# 	P_prime = Rfast::matrnorm(n, ncol(ch)) %*% ch / sqrt(nu)

# 	# V_i = crossprod(P_prime) 
# 	lambda = corpcor::estimate.lambda(P_prime, verbose=FALSE)
# 	# lambda = 0
# 	C = crossprod(scale(P_prime)) / (nrow(P_prime) - 1)
# 	V_i = (1-lambda)*C + lambda*diag(1, ncol(P_prime))
# 	V_i = diag(stders) %*% V_i %*% diag(stders)

# 	# beta_i = rmvnorm( 1, rep(0, ncol(V_i)), V_i)
# 	# beta_i = c(beta_i)

# 	ones = matrix(1,1,length(beta_i))

# 	a = solve(V_i, beta_i)
# 	b = solve(V_i, t(ones))

# 	newx <- (ones %*% a)/(ones %*% b)
#     newv <- 1/(ones %*% b)
# 	stat = (newx / sqrt(newv))^2
# 	stat[1]
# }

# # WHERE to apply lambda?
# # shrink after wishart sampling

# library(EnvStats)
# with(df, LS.emp.resid(beta, se, P, nu[1], n.mc.samples=1e3))

# with(df, remaCor::LS.empirical(beta, se, S, nu[1], n.mc.samples=1e3))


# beta = df$beta
# stders = df$se

# # Use eigen-axes
# get_sample = function(C_chol, nu){

# 	# RSS ~ W(n,Sigma)

# 	# C_chol = chol(C)
# 	# C_i = rwishart(nu, C) / nu
# 	# C_i_chol = rwishartc(nu, C) / sqrt(nu)
# 	C_i_chol = rwishart_chol(nu, C_chol) / sqrt(nu)
# 	C_i_chol = msqrt(crossprod(C_i_chol))
# 	V_i_chol = C_i_chol / sqrt(nu)

# 	beta_i = V_i_chol %*% rnorm(nrow(C_chol))

# 	C_i_chol = rwishart_chol(nu, C_chol) / sqrt(nu)
# 	# C_i_chol = msqrt(crossprod(C_i_chol))
# 	V_i_chol = C_i_chol / sqrt(nu)

# 	ones = matrix(1,1,length(beta_i))
# 	V_i = crossprod(V_i_chol)

# 	a = solve(V_i, beta_i)
# 	b = solve(V_i, t(ones))

# 	newx <- (ones %*% a)/(ones %*% b)
#     newv <- 1/(ones %*% b)
# 	(newx / sqrt(newv))^2
# }



