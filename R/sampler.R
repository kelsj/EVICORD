#' function: sampler.r
#' gibbs sampling of hierarchical model from allele pair recombination distances
#' user provided values: allele count, alphaRec, betaRec, ..., outfile

sampler = function(dists,ac,alphaRec,betaRec,alphaIBD,alpha0,beta0,piPriors,niter,outfile="gibbs_sampler_out.rds") {
	
	#allele count
	npairs = choose(ac,2)
	
	#set value for alpha IBD, intial value for beta IBD
	betaIBD_init = 0.1
	
	alphas = c(alphaIBD, alphaRec)
	betas = c(betaIBD_init, betaRec)
	
	#rescale piPriors to sum to 1
	piPriors = piPriors/sum(piPriors)
	
	#partitions from getPartitions
	part = getPartitions(ac)
	nrec_k = part$nrec_k
	nibd1_k = part$nibd1_k
	nibd2_k = part$nibd2_k
	partitions = part$partitions
	nk = part$nk
	
	#initialize beta, t, z, v, pi values
	beta_i = matrix(ncol=2,nrow=niter+1)
	beta_i[1,] = beta_truth
	t1_i = matrix(nrow=niter+1,ncol=n) #IBD grp1: rows = iters, cols = values
	t2_i = matrix(nrow=niter+1,ncol=n) #IBD grp2
	t3_i = matrix(nrow=niter+1,ncol=n) #rec pairs
	t1_i[1,] = rgamma(n,shape=alphas[1],rate=betas[1])
	t2_i[1,] = rgamma(n,shape=alphas[1],rate=betas[1])
	t3_i[1,] = rgamma(n,shape=alphas[2],rate=betas[2])
	
	#initialize z values, sample with piPriors probs
	z_chain = matrix(nrow=niter+1,ncol=n)
	z_chain[1,] = sample(1:nk,n,replace=T,prob= piPriors)
	
	#initial v values, sample from piPriors
	v_chain = matrix(nrow=niter+1,ncol=m)
	v_chain[1,] = sample(1:2,m,replace=T,prob=c(sum((npairs-nrec_k)*piPriors),sum(nrec_k*piPriors)))
	
	#initial pi values
	pi_chain = matrix(nrow=niter+1,ncol=nk)
	pi_chain[1,] = piPriors
	
	
	for (i in 1:niter){
		#sample beta for ibd
		beta_i[i+1,1] = sample_beta(c(t1_i[i,],t2_i[i,]) %>% na.omit(), alphas[1],alpha0,beta0)
		beta_i[i+1,2] = betaRec
	
		#update z,v,t
		zvt = sample_t_z_v(alphas,beta_i[i+1,],pi_chain[i,],dists$distL_cM,dists$distR_cM)
		z_chain[i+1,] = zvt$z
		v_chain[i+1,] = zvt$v
		t1_i[i+1,] = zvt$t[,1]
		t2_i[i+1,] = zvt$t[,2]
		t3_i[i+1,] = zvt$t[,3]
	
		#update pi
		pi_chain[i+1,] = sort(sample_pi(z_chain[i+1,],piPriors),decreasing=T)
	}
	
	outlist = list("nVars"=n,"beta"=beta_i,"pi"=pi_chain,"t1"=t1_i,"t2"=t2_i,"t3"=t3_i,"z"=z_chain,"v"=v_chain)
	saveRDS(outlist,file=outfile)
	
	
	return(paste("Gibbs output saved to ",outfile,sep=""))

}