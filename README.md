# EVICORD

EVICORD (EValuating IBD Consistency via Obligate Recombination Distance)

This R package implements a Gibbs sampler for the Bayesian hierarchical model described in [Johnson & Voight 2020](https://www.biorxiv.org/content/10.1101/2020.05.26.117358v1).

The purpose of this method is to classify variants as likely IBD or non-IBD using the pairwise recombination distances between allele pairs.


## Installation

```R
#install from github
library(devtools)
github_install("kelsj/ibdibsR")
#load package
library(ibdibsR)
```


## Steps to run sampler


1. Load pairwise recombination distances in centimorgans (file format should have 3 columns: varID, distL, distR). The distances should be ordered by variant and by allele pairs. Within each variant, the order fo allele pair recombination distances matters; e.g. for allele count 3, the pairwise distances should be in the order 1-2,1-3,2-3; for allele count 4, the order should be 1-2,1-3,1-4,2-3,2-4,3-4; etc.

```R
dists = ibdibsR::load_dists("sim_dists_ex.txt",header=T,ac=4)
```

2. Run the Gibbs sampler on the pairwise recombination distances. Several parameter values must be provided at this step: the allele count (**ac**), values for alpha and beta for recurrent/non-IBD variants (**alphaRec**, **betaRec**; can be calculated with nonIBD\_param\_values.R, see below), alpha for IBD variants (**alphaIBD**; recommended to try multiple values and see how this affects your results), priors for beta for IBD variants (**alpha0**, **beta0**), vector of priors for pi (**piPriors**, should have as many values as there are categories of k; i.e. 1+number of non-IBD partitions), and the number of iterations to run (**niter**).

```R
ibdibsR::sampler(dists,ac=4,alphaRec=40,betaRec=0.1,alphaIBD=10,alpha0=1,beta0=10,piPriors=c(0.4,0.4,0.2),niter=100,outfile="gibbs_sampler_out.rds") #returns RDS object of output: gibbs_sampler_out.rds
```

3. Visualize results from gibbs_sampler_out.rds. Required parameters: **filename**, e.g. gibbs\_sampler\_out.rds; **ac**, allele count; **tf** integer to thin chains by for visualization. Output plot panels:

    1. chain of beta samples by iteration
    2. chain of pi samples by iteration
    3. density of t estimates for each allele pair category (t1,t2,t3)
    4. empirical cumulative density of the fraction of z samples for each variant greater than x (with x=1,2,3,...); i.e., zGT1 = the fraction of posterior samples that the variant is assigned to be non-IBD.

```R
ibdibsR::viz_res(filename="gibbs_sampler_out.rds",ac=4,tf=5,outfile="gibbs_sampler_out_viz.pdf") #saves visualizations in pdf gibbs_sampler_out_viz.pdf
```

## Estimating beta for non-IBD variants from data

Posterior values of alpha and beta for the distribution of TMRCAs for non-IBD variants need to be supplied to run the sampler function above. These can be estimated from the pairwise recombination distances of non-IBD allele pairs from multiallelic or known non-IBD variants using the function nonIBD\_param\_values. Requires input of non-IBD variants' pairwise recombination distance (columns: varID, distL, distR). Several parameter values must be provided: the number of non-IBD allele pair distances provided for each variant (**n_nonIBDpairs**), the value of alpha for non-IBD allele pairs (**alpha**); priors for beta (**alpha0**, **beta0**), the number of iterations to run (**niter**), and the output file name (**outfile**).

```R
ibdibsR::nonIBD_param_values(dists,n_nonIBDpairs=3,alpha=40,alpha0=1,beta0=10,niter=500,outfile="sampled_nonIBD_param_values_alpha40") #saves sampled values to RDS file sampled_nonIBD_param_values.rds
```



