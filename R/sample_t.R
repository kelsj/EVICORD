#' function: sample t (TMRCA)
#' sample t from gamma distribution, given distance d, alpha, beta
#' t ~ gamma(alpha,beta)

sample_t = function(d,alpha,beta){
	#d = vector of distances
	t = rgamma(1,shape=(alpha+length(d)),rate=(beta+sum(d)/50))
	return(t)
}
