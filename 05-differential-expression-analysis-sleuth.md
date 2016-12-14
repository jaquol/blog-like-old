## 05. Differential expression analysis (DEA) with Sleuth

I use [Kallisto](https://pachterlab.github.io/kallisto/about.html)'s transcript abundances to perform DEA using [Sleuth](http://pachterlab.github.io/sleuth/).

In the typical Sleuth's workflow, the likelihood ratio test (LRT) is applied. Briefly, the LRT models the likelihood of the data given 2 models:
- *full*: transcript abundance affected on one or more dependent variables (here just being treated or not)
- *reduced*: transcript abundance unaffected by the treatment (null hypothesis)

The LRT then estimates for each transcript the ratio of the 2 likelihoods and produces a q-value (i.e. p-value adjusted for multiple testing by means of false discovery rate, FDR) which can be used as measure of significance. Note that the LRT does not produce any metric equivalent to the fold change, which indicates the maginitude of the change in expression between the 2 conditions and is commonly reported in DEA. That's why the script below also applies the Wald test (WT), which is somewhat related to the LRT and is also used to test for differential expression. WT is used becase it generates the beta statistic, which approximates to the log2 fold change in expression between the 2 condition tested. However, LRT is considerd a better test than the WT (see [here](http://www.ats.ucla.edu/stat/mult_pkg/faq/general/nested_tests.htm)) and thus significance filtering is based on LRT's q-values.

My implementation of Sleuth (see below) is largerly inspired on the [introduction to Sleuth](https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html) and on this [blog post](http://achri.blogspot.com.es/2016/10/using-kallisto-sleuth.html) explaining the usage of the WT.

When applying the WT, in order to get beta, **understanding how Sleuth reads the 2 treatments/conditions tested is key to understand how expression changes between them**.

Imagine that we have 3 samples for each of 2 conditions (untreated and treated):

| sample_id  | condition  | path  |
|---|---|---|
| sample01  | untreated  | path1  |
| sample02  | untreated  | path2  |
| sample03  | untreated  | path3  |
| sample04  | treated  | path4  |
| sample05  | treated  | path5  |
| sample06  | treated  | path6  |

At some point we load the Kallisto-processed data and make a regression model using 'condition' as the dependent variable:
```
so <- sleuth_prep(tab_metadata, ~ condition)
```
And then we apply the WT:
```
so <- sleuth_wt(so, paste0('condition'))
```
Internally, with `sleuth_prep` Sleuth will transform elements in the condition field to 0s and 1s **in alphabetical order** and then WT's beta values will be relative to the 0 condition; that is, positive beta values showing transcripts in which expression is greater in condition 1 than in condition 0. But this would result in counter-intuitive results in our experiment because beta values would be relative to the '**T**reated' samples. So the trick is making that the *reference* condition ranks first alphabetically speaking. Below is my solution to avoid having to play excessively with the condition names.

**The code** 

```
diff_exp_sleuth <- function(condition1, condition2) {

	# re-format metadata to meet sleuth's required format
	# later on, the 'sleuth_prep' function of sleuth will assign 0=untreated and 1=treated condition
	# importantly, by default such function will set '0' to first alphabetical match in the values
	# found in the 'condition' field so here we force that 0=condition1
	tab_metadata <- metadata
	rownames(tab_metadata) <- tab_metadata$sample_id
	
	conditions <- c()
	for (i in 1:nrow(tab_metadata)) {
		tag <- paste(tab_metadata[i, 'treatment'], tab_metadata[i, 'treatment_time'], sep = '_')
		if (tag == condition1) {
			conditions <- c(conditions, 0)
		}
		else if (tag == condition2) {
			conditions <- c(conditions, 1)
		}
		else {
			conditions <- c(conditions, -1)
		}
	}
 	tab_metadata$condition <- conditions

	# subset samples to include only those from the 2 conditions compared
	cond1 <- tab_metadata$condition == 0
  	cond2 <- tab_metadata$condition == 1
 	tab_metadata <- tab_metadata[cond1 | cond2, ]

 	samples <- tab_metadata$sample_id

 	# get paths to data
	paths <- c()
	for (s in samples) {
		p <- paste0(SAMPLES, "/", s, "/quantifications/kallisto/", assembly_version, '/', sequencing_type)
		paths <- c(paths, p)
	}
	tab_metadata$path <- paths

 	# rename to meet sleuth's format requirements
 	tab_metadata <- tab_metadata[, c('sample_id', 'condition', 'path')]
 	names(tab_metadata)[1] <- 'sample'

	 # (1) load the kallisto processed data and make a regression model using 'condition' as the dependent variable
  	so <- sleuth_prep(tab_metadata, ~ condition)
  
	# (2) estimate parameters for the sleuth response error measurement (full) model as responding to the 'condition' factor
  	so <- sleuth_fit(so)
  
	# (3) Create another model where the gene expression is not dependent on any factor.
  	so <- sleuth_fit(so, ~1, 'reduced')
  
	# (4.1) Run a likelihood ratio test (LRT) between the two models to see what transcripts appear 
	# to really be affected by the time factor value
	so <- sleuth_lrt(so, 'reduced', 'full')

	# (4.2) Run the Wald test (WT), a statistical tests which:
	# - is somewhat related to the LRT and is also used to test for differential expression
	# - LRT is considerd a better test than the WT but
	# - WT is used becase it generates the beta statistic, which approximates to the fold change in expression between
	# the 2 condition tested, which is typically reported in differential expression analysis
	so <- sleuth_wt(so, paste0('condition'))

	# export normalised abundance values
	condition1_name <- tolower(condition1)
	condition2_name <- tolower(condition2)
	otab = paste0(ANALYSIS, "/tables/normalized_abundance_transcript_level_sleuth_", condition1_name, "_", condition2_name, ".tsv")
 	write.table(kallisto_table(so), otab, sep = "\t", quote = FALSE, row.names = FALSE)

	# add beta (b), beta's standard error (se_b) and the mean expression in the samples (mean_obs)
 	res_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
 	res_wt <- sleuth_results(so, 'condition')
 	res <- merge(res_lrt, res_wt[, c('target_id', 'b', 'se_b', 'mean_obs')], on = 'target_id', sort = FALSE)
	
	# export
	condition1_name <- tolower(condition1)
	condition2_name <- tolower(condition2)
	otab = paste0(ANALYSIS, "/tables/differential_expression_analysis_transcript_level_sleuth_", condition1_name, "_", condition2_name, ".tsv")
 	write.table(res, otab, sep = "\t", quote = FALSE, row.names = FALSE)

 	return(so)

}

sleuth_object = diff_exp_sleuth("Untreated_0", "Doxycycline_480")
sleuth_object = diff_exp_sleuth("Untreated_0", "Doxycycline_960")
```
