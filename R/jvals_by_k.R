#' function: jvals_by_k
#' generates matrices of all possible assignments j for each allele pair, for a partition k

#jvals for each k:
jvals_by_k = function(nk,npairs,partitions){
	jvals_k = matrix(nrow=nk,ncol=npairs) #1=IBD pair, 2=Rec pair
	grp_k = matrix(nrow=nk,ncol=npairs) #1=IBD1, 2=IBD2, 3=Rec
	for (sk in 1:nk){
		if (sk==1){
			#all pairs are IBD
			jvals_k[sk,] = rep(1,npairs)
			grp_k[sk,] = rep(1,npairs)
		} else {
			#pairs are rec/IBD depending on partition
			part = partitions[sk-1]
			jvalsRec = c()
			grpsRec = c()
			grp1 = c(1:part)
			grp2 = c((part+1):ac)
			for (ind1 in 1:(ac-1)){
				for (ind2 in (ind1+1):ac){
					if (ind1 %in% grp1 & ind2 %in% grp1) {
						#IBD1 pair
						jvalsRec = c(jvalsRec,1)
						grpsRec = c(grpsRec,1)
					} else if (ind1 %in% grp2 & ind2 %in% grp2) {
						#IBD2 pair
						jvalsRec = c(jvalsRec,1)
						grpsRec = c(grpsRec,2)
					} else {
						#rec pair
						jvalsRec = c(jvalsRec,2)
						grpsRec = c(grpsRec,3)
					}
				}
			}
			jvals_k[sk,] = jvalsRec
			grp_k[sk,] = grpsRec
		}
	}
	return(list("jvals_k"= jvals_k, "grp_k"= grp_k))
}
