#' getPartitions
#' from allele count (ac), generate list of possible recurrent partitions
#' nrec_k = number of recurrent allele pairs for partition k
#' nibd1_k = number of IBD allele pairs in IBD grp1 for partition k
#' nibd2_k = number of IBD allele pairs in IBD grp1 for partition k
#' partitions = partition (e.g. 1 for 1:(ac-1), 2 for 2:(ac-2))
#' nk = number of partitions for this allele count

getPartitions = function(ac) {
	npairs = choose(ac,2)
	nrec_k = c(0)
	nibd1_k = c(npairs)
	nibd2_k = c(0)
	partitions = c()
	for (part in 1:as.integer(ac/2)){
		partitions = c(partitions, part)
		nrec_k = c(nrec_k,part*(ac-part))
		nibd1_k = c(nibd1_k,choose(part,2))
		nibd2_k = c(nibd2_k,choose((ac-part),2))
	}
	nk = length(nrec_k)
	return(list("nrec_k"=nrec_k, "nibd1_k"=nibd1_k, "nibd2_k"= nibd2_k, "partitions"=partitions, "nk"=nk))
}