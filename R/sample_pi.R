#' function: sample_pi
#' sampe pi from dirichlet given assignments of k for each variant (vector z), priors (vector piPriors)
#' nk = number of k values (1=IBD, 1=1:(n-1), 2=2:(n-2), etc.)
#' @export

#' @import gtools

sample_pi = function(z,piPriors,nk){
	counts = c()
	for (c in 1:nk){
		counts = c(counts,length(z[z==c]))
	}
	counts = counts+piPriors #add priors
	return(rdirichlet(1,counts)) #pi_t: vector of length nk
}
