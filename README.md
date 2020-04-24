# ibdibsR

This R package implements a Gibbs sampler for the Bayesian hierarchical model described in [Johnson & Voight 2020](hyperlink).

The purpose of this method is to classify variants as likely IBD or non-IBD using the pairwise recombination distances between allele pairs.

The steps are as follows:

1. Calculate pairwise recombination distances from VCF formatted file and genetic map file:

```R
vcf = read.table("data.vcf",header=T)
calc_recomb_dists(vcf) #saves recomb dists in file recomb_dists.txt.gz
```

2. Load pairwise recombination distances in centimorgans (file format should have 3 columns: varID, distL, distR). If you are not using a distances file calculated from calc_recomb_dists, it is important to note that the order of the recombination distances matters; e.g. for allele count 3, the pairwise distances should be in the order 1-2,1-3,2-3.

```R
dists = load_dists_file("recomb_dists.txt.gz", head=T, ac=3)
```

3. 