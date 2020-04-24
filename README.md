# ibdibsR

This R package implements a Gibbs sampler for the Bayesian hierarchical model described in [Johnson & Voight 2020](hyperlink).

The purpose of this method is to classify variants as likely IBD or non-IBD using the pairwise recombination distances between allele pairs.


## Steps to run sampler

The steps are as follows:

1. Calculate pairwise recombination distances from VCF formatted file and genetic map file:

```R
vcf = read.table("data.vcf",header=T)
map = read.table("map.txt",header=T)
calc_recomb_dists(vcf,map) #saves recomb dists in file recomb_dists.txt.gz
```

2. Load pairwise recombination distances in centimorgans (file format should have 3 columns: varID, distL, distR). If you are not using a distances file calculated from calc\_recomb\_dists, it is important to note that the order of the recombination distances matters; e.g. for allele count 3, the pairwise distances should be in the order 1-2,1-3,2-3; for allele count 4, the order should be 1-2,1-3,1-4,2-3,2-4,3-4; etc.

```R
dists = load_dists_file("recomb_dists.txt.gz", head=T, ac=3)
```

3. Run the Gibbs sampler on the pairwise recombination distances. Several parameter values must be provided at this step: the allele count (**ac**), values for alpha and beta for recurrent/non-IBD variants (**alphaRec**, **betaRec**; can be calculated with nonIBD\_param\_values.R, see below), alpha for IBD variants (**alphaIBD**; recommended to try multiple values and see how this affects your results), priors for beta for IBD variants (**alpha0**, **beta0**), vector of priors for pi (**piPriors**, should have as many values as there are categories of k; i.e. 1+number of non-IBD partitions), and the number of iterations to run (**niter**).

```R
sampler(dists,ac,alphaRec,betaRec,alphaIBD,alpha0,beta0,piPriors,niter) #returns RDS object of output: gibbs_sampler_out.rds
```

4. Visualize results from gibbs_sampler_out.rds. Required parameters: **filename**, e.g. gibbs\_sampler\_out.rds; **ac**, allele count; **tf** integer to thin chains by for visualization. Output plot panels:

  1. chain of beta samples by iteration
  2. chain of pi samples by iteration
  3. density of t estimates for each allele pair category (t1,t2,t3)
  4. empirical cumulative density of the fraction of z samples for each variant greater than x (with x=1,2,3,...); i.e., zGT1 = the fraction of posterior samples that the variant is assigned to be non-IBD.

```R
viz_res(filename="gibbs_sampler_out.rds") #saves visualizations in pdf: gibbs_sampler_out_viz.pdf
```

## Estimating alpha and beta for non-IBD variants from data

Posterior values of alpha and beta for the distribution of TMRCAs for non-IBD variants need to be supplied to run the sampler function above. These can be estimated from the pairwise recombination distances of non-IBD allele pairs from multiallelic or known non-IBD variants using the function nonIBD\_param\_values. 





