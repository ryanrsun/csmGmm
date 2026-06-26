# csmGmm Description and quick-start

This R package implements the conditionally symmetric multidimensional Gaussian mixture model (csmGmm) for large-scale testing of composite null hypotheses in genetic association applications such as replication analysis, pleiotropy analysis, and mediation analysis. In such analyses, we typically have K sets of test statistics corresponding to the same J SNPs, where K is a small number (e.g. 2 or 3) and J is large (e.g. 1 million). For each SNP, we want to know if we can reject all K individual nulls. The methodology paper reference is: Sun R, McCaw ZR, Lin X. Testing a large number of composite null hypotheses using conditionally symmetric multidimensional gaussian mixtures in genome-wide studies. Journal of the American Statistical Association. 2025 Apr 3;120(550):605-17.

For a quick-start replication analysis, try the following:

```r
# Load library
library(csmGmm)

# Download and read data
dat <- read_gwas_set(c(
  "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010098/GCST010098.tsv",
  "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90551001-GCST90552000/GCST90551892/GCST90551892.tsv.gz"
))

# load datasets into list
dat1 <- dat$gwas1
dat2 <- dat$gwas2
mydatasets <- list(dat1, dat2)

# merge and clean data
prepped <- prepare_csmgmm_data(datasets=mydatasets, dataset_names = NULL, z_cap = 8.1) 

# generate suggested initial parameters
initParams <- generate_init_lists(2)

# perform replication analysis
res <- symm_fit_ind_EM(testStats = prepped$clean_dat, 
                       initMuList = initParams$initMuList, 
                       initPiList = initParams$initPiList, 
                       sameDirAlt = TRUE, 
                       eps = 10^(-5), 
                       checkpoint = TRUE)

# process raw results
processed <- process_lfdr_results(orig_data=prepped$extra_dat, 
                                  lfdrResults=res$lfdrResults, fdr_threshold = 0.1) 

# produce visualizations
figs <- create_plots(plotData = processed$all_snps_info)
```
