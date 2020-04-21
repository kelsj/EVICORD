library(dplyr)
library(ggplot2)
library(invgamma)
library(gtools)
library(gridExtra)
library(cowplot)
library(tidyr)
library(stringr)

#REP2 in output file

#gibbs sample to estimate hierarchical model including t estimates for each var
#use model for recurrent pairs alpha/beta previously estimated from rec dists only:
# alphaRec = 140
# betaRec = 0.328
alphaRec = 40
betaRec = 0.0859
# alphaRec = 30
# betaRec = 0.0624
# alphaRec = 20
# betaRec = 0.0392
# alphaRec = 10
# betaRec = 0.0170
niter = 5000

args = commandArgs(TRUE)
if (length(args)==0){
	stop("input: Rscript hier_model_fitT_uk10k.R filename AC")
}

filename = args[1]
ac = as.integer(args[2])
npairs = choose(ac,2)

filesub = str_split(filename,"recDists.txt.gz")[[1]][1]

#read inds that pass QC
passInds = read.table("/project/voight_subrate/kelsj/uk10k/recDists/alspac_twins/classify_by_dists_v4/probFiles/ALSPAC_TWINS_3621samples_PASS_UK10K_PROJECT_QC_PUBLICIDS.txt",header=F) %>% tbl_df()

#read input
#file format:
#varID ind1 ind2 distL_bp distR_bp distL_cM distR_cM endL endR
dists = read.table(gzfile(paste("/project/voight_subrate/kelsj/uk10k/recDists/alspac_twins/parallelize/",filename,sep="")),header=T) %>% tbl_df()

#filter dists for those inds that pass UK10K QC
#filter recomb distances of zero (@ ends of chroms)
dists = dists %>% filter(ind1 %in% passInds$V1 & ind2 %in% passInds$V1, distL_cM>0, distR_cM>0)

#filter for those with correct # allele pairs
incl_varIDs = dists %>% group_by(varID) %>% summarise(count=n()) %>% filter(count == npairs)
dists = dists %>% filter(varID %in% incl_varIDs$varID)

#number of variants included
n = length(unique(dists$varID))
m = n*npairs

#k labels: k=1 IBD, k=2 rec1:n-1, k=3 rec2:n-2, ...
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

#given ac, list of inds in one cluster, assign ibd/rec pairs
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
#matrices of ibd pairs for each rec k
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

#jvals for each k:
jvals_k = matrix(nrow=nk,ncol=npairs) #1=IBD pair, 2=Rec pair
grp_k = matrix(nrow=nk,ncol=npairs) #1=IBD1, 2=IBD2, 3=Rec
for (sk in 1:nk){
	if (sk==1){
		#all IBD
		jvals_k[sk,] = rep(1,npairs)
		grp_k[sk,] = rep(1,npairs)
	} else {
		#rec/IBD depending on partition
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

#sample t ~ gamma(alpha,beta)
sample_t = function(d,alpha,beta){
	#d = vector of distances
	t = rgamma(1,shape=(alpha+length(d)),rate=(beta+sum(d)/50))
	return(t)
}

#sample beta from beta|d,t,alpha
sample_beta = function(t,alpha,alpha0,beta0){
	#t = vector of ages
	beta = rgamma(1,shape=(alpha0+length(t)*alpha),rate=(beta0+sum(t)))
	return(beta)
}

#z: assignments for k for each allele pair, vector length n
sample_pi = function(z,piPriors){
	counts = c()
	for (c in 1:nk){
		counts = c(counts,length(z[z==c]))
	}
	counts = counts+piPriors #add priors
	return(rdirichlet(1,counts)) #pi_t: vector of length nk
}

#z: assignments for k for each allele pair, vector length n
#alpha = c(alphaIBD,alphaRec), beta=c(betaIBD,betaRec), pi=c(p1,p2,..,pk), x_i1=distsL, x_i2=distsR
sample_t_z_v = function(alpha,beta,pi,x_i1,x_i2){
#	zvt = sample_t_z_v(alpha_i[i+1,],beta_i[i+1,],pi_chain[i+1,],dists$distL_cM,dists$distR_cM)
	#every npairs are from one variant
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
			#all pairs j=1 (IBD)
			sj = rep(1,npairs)
			v[(d*npairs-(npairs-1)):(d*npairs)] = sj
			
			## save t IBD ##
			t[d,] = c(tIBD,NA,NA)
			
		} else {
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

#"true" parameters for alpha, beta, pi:
alpha_truth = c(10, alphaRec)
beta_truth = c(0.1,betaRec)
if (ac>=4){
	pi_truth = c(0.75,0.15,rep(0.1/(nk-2),(nk-2))) #true proportions of k components: 75% IBD, 10% 1:(n-1), 15% the rest
} else {
	pi_truth = c(0.75,0.25) #true proportions of k components (ac 2 or 3, nk==2)
}

#gibbs sampling for alpha, beta, t, & k from d observations
beta_i = matrix(ncol=2,nrow=niter+1)
beta_i[1,] = beta_truth
t1_i = matrix(nrow=niter+1,ncol=n) #IBD grp1: rows = iters, cols = values
t2_i = matrix(nrow=niter+1,ncol=n) #IBD grp2
t3_i = matrix(nrow=niter+1,ncol=n) #rec pairs
t1_i[1,] = rgamma(n,shape=alpha_truth[1],rate=beta_truth[1])
t2_i[1,] = rgamma(n,shape=alpha_truth[1],rate =beta_truth[1])
t3_i[1,] = rgamma(n,shape=alpha_truth[2],rate =beta_truth[2])
z_chain = matrix(nrow=niter+1,ncol=n)
z_chain[1,] = sample(1:nk,n,replace=T,prob=pi_truth)
v_chain = matrix(nrow=niter+1,ncol=m)
v_chain[1,] = sample(1:2,m,replace=T,prob=c(sum((npairs-nrec_k)*pi_truth),sum(nrec_k*pi_truth)))
pi_chain = matrix(nrow=niter+1,ncol=nk)
pi_chain[1,] = pi_truth

#pi priors
piPriors = pi_truth*10
#beta priors
alpha0 = 1
beta0 = 10

for (i in 1:niter){
	#print(i)
	#sample beta for ibd
	beta_i[i+1,1] = sample_beta(c(t1_i[i,],t2_i[i,]) %>% na.omit(), alpha_truth[1],alpha0,beta0)
	beta_i[i+1,2] = betaRec

	#update z,v,t
	zvt = sample_t_z_v(alpha_truth,beta_i[i+1,],pi_chain[i,],dists$distL_cM,dists$distR_cM)
	z_chain[i+1,] = zvt$z
	v_chain[i+1,] = zvt$v
	t1_i[i+1,] = zvt$t[,1]
	t2_i[i+1,] = zvt$t[,2]
	t3_i[i+1,] = zvt$t[,3]

	#update pi
	pi_chain[i+1,] = sort(sample_pi(z_chain[i+1,],piPriors),decreasing=T)
}

outlist = list("nVars"=n,"beta"=beta_i,"pi"=pi_chain,"t1"=t1_i,"t2"=t2_i,"t3"=t3_i,"z"=z_chain,"v"=v_chain)
saveRDS(outlist,file=paste(filesub,"hierModel_estT_gammaPrior_alphaRec",alphaRec,"_rep2.rds",sep=""))



