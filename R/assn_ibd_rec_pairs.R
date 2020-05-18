#' assn_ibd_rec_pairs
#' given ac and vector of inds in a cluster (cluster1), get inds in each ibd/rec pair
#' @export

assn_ibd_rec_pairs = function(ac,cluster1){
	cluster2 = c(1:ac)[!(c(1:ac) %in% cluster1)]
	#IBD pairs those in same cluster, rec across clusters
	ibd1Pairs = c()
	ibd2Pairs = c()
	recPairs = c()
	pairNum = 0
	for (i1 in 1:(ac-1)){
		for (i2 in (i1+1):ac){
			pairNum = pairNum+1
			if ((i1 %in% cluster1) & (i2 %in% cluster1)){
				ibd1Pairs = c(ibd1Pairs,pairNum)
			} else if ((i1 %in% cluster2) & (i2 %in% cluster2)) {
				ibd2Pairs = c(ibd2Pairs,pairNum)
			} else {
				recPairs = c(recPairs,pairNum)
			}
		}
	}
	return(list("ibd1p"=ibd1Pairs,"ibd2p"=ibd2Pairs,"recp"=recPairs))
}
