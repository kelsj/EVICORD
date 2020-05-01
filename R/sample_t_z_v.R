#' function: sample_t_z_v
#' z = assignments of k for each variant
#' t = assignments of t for each variant
#' v = assignments of j for each allele pair

#' input:
#' alpha = c(alphaIBD,alphaRec)
#' beta=c(betaIBD,betaRec)
#' pi=c(p1,p2,..,pk)
#' x_i1=distsL
#' x_i2=distsR
#' distances are in order so each npairs are for one variant
#' npairs = numbr of allele pairs for these allele count (=choose(ac,2))

sample_t_z_v = function(alpha,beta,pi,x_i1,x_i2,ac){
	npairs = choose(ac,2)
	
	#get matrices of pairs for each pair type
	pairmat = ibdRecPairMat(ac)
	ibd1AssByK = pairmat$ibd1AssByK
	ibd2AssByK = pairmat$ibd2AssByK
	recAssByK = pairmat$recAssByK
	
	z = c()
	v = c()
	t = matrix(nrow=n,ncol=3)
	for (d in 1:n){
		#arrange by ind pairs 1-2,1-3,...,1-n,2-3,2-4,etc.
		xi1 = x_i1[(d*npairs-(npairs-1)):(d*npairs)]
		xi2 = x_i2[(d*npairs-(npairs-1)):(d*npairs)]
		
		#estimate tIBD for each allele pair, p(k==1)
		tIBD = sample_t(c(xi1,xi2),alpha=alpha[1], beta=beta[1])
		#start vector of pk's with p(k|d,t,alpha,beta,pi)
		pk = c(prod(dexp(c(xi1,xi2),rate=tIBD/50))*dgamma(tIBD,shape=alpha[1],rate=beta[1])*pi[1])
		#if pk_ibd==0, v long distances, call as IBD
		if (pk[1]==0){
			z[d] = 1
			t[d,] = c(tIBD,NA,NA)
		} else {		
			#iterate rec k vals
			pk_j_poss = list()
			t_j_poss = list()
			for (part in 1:as.integer(ac/2)){
				#go through possible rec/ibd assignments
				ibd1Ass = ibd1AssByK[[paste(part+1)]]
				ibd2Ass = ibd2AssByK[[paste(part+1)]]
				recAss = recAssByK[[paste(part+1)]]
				#store pk, t samples for each peralphatation
				pk_poss = c()
				t_poss = matrix(nrow=nrow(recAss),ncol=3)
				for (q in 1:nrow(recAss)){
					#sample t for each cluster
					tIBD1 = ifelse(ncol(ibd1Ass)>0,sample_t(c(xi1[ibd1Ass[q,]],xi2[ibd1Ass[q,]]),alpha=alpha[1],beta=beta[1]),NA)
					tIBD2 = ifelse(ncol(ibd2Ass)>0,sample_t(c(xi1[ibd2Ass[q,]],xi2[ibd2Ass[q,]]),alpha=alpha[1],beta=beta[1]),NA)
					tRec = sample_t(c(xi1[recAss[q,]],xi2[recAss[q,]]),alpha=alpha[2],beta=beta[2])
					#p(d|t)
					pIBD1d = ifelse(ncol(ibd1Ass)>0,prod(dexp(c(xi1[ibd1Ass[q,]],xi2[ibd1Ass[q,]]),rate=tIBD1/50)),1)
					pIBD2d = ifelse(ncol(ibd2Ass)>0,prod(dexp(c(xi1[ibd2Ass[q,]],xi2[ibd2Ass[q,]]),rate=tIBD2/50)),1)
					pRecd = prod(dexp(c(xi1[recAss[q,]],xi2[recAss[q,]]),rate=tRec/50))
					#p(t|k,alpha,beta)
					pIBD1t = dgamma(tIBD1,shape=alpha[1],rate=beta[1])
					pIBD2t = dgamma(tIBD2,shape=alpha[1],rate=beta[1])
					pRect = dgamma(tRec,shape=alpha[2],rate=beta[2])
					pt = mean(c(pIBD1t,pIBD2t,pRect),na.rm=T)
					pk_poss[q] = pIBD1d*pIBD2d*pRecd*pt
					t_poss[q,] = c(tIBD1,tIBD2,tRec)
				}
				pk[(part+1)] = mean(pk_poss*pi[part+1]) #weight pk_poss by pi for this k
				#store j assignment probs
				pk_j_poss[[paste(part+1)]] = pk_poss
				#store j t samples
				t_j_poss[[paste(part+1)]] = t_poss
			}
			z[d] = sample(1:nk,1,prob=pk,replace=T)
		}
		if (z[d]==1){
			#IBD variant, all pairs j=1
			sj = rep(1,npairs)
			v[(d*npairs-(npairs-1)):(d*npairs)] = sj
			
			## save t IBD ##
			t[d,] = c(tIBD,NA,NA)
			
		} else {
			#recurrent variant, use possible j assignments
			ncombs = nrow(recAssByK[[paste(z[d])]])
			#sample a combination
			#divide probs by pi for chosen k b/c prob of k calculated above was weighted by pi_t
			sc = sample(1:ncombs,1,prob=pk_j_poss[[paste(z[d])]],replace=T)
			#create ordered vector of j assignments
			sj = rep(0,npairs)
			sj[ibd1AssByK[[paste(z[d])]][sc,]] = 1 #IBD pairs
			sj[ibd2AssByK[[paste(z[d])]][sc,]] = 1 #IBD pairs
			sj[recAssByK[[paste(z[d])]][sc,]] = 2 #rec pairs
			v[(d*npairs-(npairs-1)):(d*npairs)] = sj
			
			## save t Rec samples ##
			t[d,] = t_j_poss[[paste(z[d])]][sc,]						
		}
	}
	return(list("z"=z,"v"=v,"t"=t))
}

