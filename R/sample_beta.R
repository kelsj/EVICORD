#' function: sample_beta
#' sample beta from beta|d,t,alpha

sample_beta = function(t,alpha,alpha0,beta0){
	#t = vector of ages
	beta = rgamma(1,shape=(alpha0+length(t)*alpha),rate=(beta0+sum(t)))
	return(beta)
}
