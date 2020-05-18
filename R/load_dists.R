#' function: load_dists
#' user provides recombination distances file, if file contains header(T/F), allele count (ac)
#' distances file format: varID distL_cM distR_cM
#' this script checks that there are 3 cols, calculate # allele pairs
#' filter for variants with correct # of allele pairs
#' @export

load_dists = function(distsFile,header=T,ac){
	dists = read.table(distsFile,header=header) %>% tbl_df()
	
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

	return(dists)
}

