#' function: nonIBD_param_values
#' gibbs sampling of t, beta from recomb distances of known non-IBD allele pairs 
#' input: dists in format varID, distL, distR
#' user provided values: n_nonIBDpairs, alpha, alpha0, beta0, niter

#sample t from t|d,alpha,beta
sample_t_nonIBD = function(d,alpha,beta){
	#d = vector of distances
	n = length(d)
	t = rgamma(1,shape=(alpha+n),rate=(beta+sum(d)/50))
	return(t)
}

#sample beta from beta|d,t,alpha
sample_beta_nonIBD = function(t,alpha,alpha0,beta0){
	#t = vector of ages
	n = length(t)
	beta = rgamma(1,shape=(alpha0+n*alpha),rate=(beta0+sum(t)))
	return(beta)
}

#gibbs sampling for alpha, beta, & t from d observations
nonIBD_param_values = function(dists,n_nonIBDpairs,alpha,alpha0,beta0,niter,outfile="sampled_nonIBD_param_values") {	
	#filter to correct # non-IBD pairs per variant
	varsKeep = dists %>% group_by(varID) %>% summarise(n=n()) %>% filter(n==n_nonIBDpairs)
	dists = dists %>% filter(varID %in% varsKeep$varID)
	nvar = nrow(varsKeep)
	
	#reshape dists: cols=variants, rows=pairwise dists
	dL = matrix(dists[,2],ncol=nvar,nrow=n_nonIBDpairs)
	dR = matrix(dists[,3],ncol=nvar,nrow=n_nonIBDpairs)
	dists = rbind(dL,dR)
	
	#initial param values
	beta_i = rep(0,niter+1)
	beta_i[1] = 0.1
	t_i = matrix(nrow=niter,ncol=nvar)
	
	#sample
	for (i in 1:niter){
		#sample t for each variant
		t_i[i,] = apply(dists,2,sample_t_nonIBD,alpha,beta=beta_i[i])
		#sample beta
		beta_i[i+1] = sample_beta_nonIBD(t=t_i[i,], alpha,alpha0,beta0)
	}
	#save output
	outdat = list("betaSamp"=beta_i, "tSamp"=t_i)
	saveRDS(outdat,file=paste(outfile,".rds",sep=""))
	return(paste("output saved in ",outfile,".rds",sep=""))
	
}

