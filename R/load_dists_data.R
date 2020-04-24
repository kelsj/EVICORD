#' function: load_dists_data
#' user provides allele count (ac) and recombination distance data in format:
#' varID distL_cM distR_cM
#' this script checks that there are 3 cols, calculate # allele pairs
#' filter for variants with correct # of allele pairs

load_dists = function(dists,ac){
	#check dimensions of data
	ncol = dim(dists)[2]
	if (ncol!=3){
		return("input error: need 3 cols (varID, distL_cM, distR_cM)")
	}
	
	#filter for those with correct # allele pairs
	npairs = choose(ac,2)
	colnames(dists) = c("varID", "distL_cM", "distR_cM")
	incl_varIDs = dists %>% group_by(varID) %>% summarise(count=n()) %>% filter(count == npairs)
	dists = dists %>% filter(varID %in% incl_varIDs$varID)

	#number of variants included
	n = length(unique(dists$varID))
	m = n*npairs

	return(dists,n,m)
}

