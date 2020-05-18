#' load recombination distances file
#'
#' reads in distances file and makes sure it is in the correct format (columns: varID distL_cM distR_cM)
#' checks that there are 3 cols, calculate # allele pairs
#' filter for variants with correct # of allele pairs
#'
#' @param distsFile recombination distances file
#' @param header logical, does file contains header(T/F)
#' @param ac allele count
#'
#' @return tbl of filtered recombination distances
#' 
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

