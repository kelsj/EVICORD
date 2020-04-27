#' ibdRecPairMat
#' for each recurrent partition k for a given allele count (ac), get matrices of all possible IBD/rec pairs

ibdRecPairMat = function(ac){
	ibd1AssByK = list()
	ibd2AssByK = list()
	recAssByK = list()
	for (part in 1:as.integer(ac/2)){
		j_comb_k = combn(ac,part)
		ibd1Ass = matrix(nrow=ncol(j_comb_k),ncol=nibd1_k[(part+1)])
		ibd2Ass = matrix(nrow=ncol(j_comb_k),ncol=nibd2_k[(part+1)])
		recAss = matrix(nrow=ncol(j_comb_k),ncol=nrec_k[(part+1)])
		for (q in 1:ncol(j_comb_k)){
			pair_assign = assn_ibd_rec_pairs(ac,j_comb_k[,q])
			ibd1Ass[q,] = pair_assign$ibd1p
			ibd2Ass[q,] = pair_assign$ibd2p
			recAss[q,] = pair_assign$recp
		}
		ibd1AssByK[[paste(part+1)]] = ibd1Ass
		ibd2AssByK[[paste(part+1)]] = ibd2Ass
		recAssByK[[paste(part+1)]] = recAss
	}
	return(list("ibd1AssByK"= ibd1AssByK, "ibd2AssByK"= ibd2AssByK, "recAssByK"= recAssByK))
}
