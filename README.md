# DesignTwoProportions

This repository contains R code associated with the article titled "Are the Tests Overpowered or Underpowered? A Unified Solution to Correctly Specify Type I Errors in the Design of Clinical Trials for Two Sample Proportions." This article was published in *Statistics in Medicine*. Please fell free to email peiran.liu@uconn.edu if you have any questions or comments.

The codes aid in the design of two-sample proportion tests using precise calculations. The repository includes two folders: one contains the functions for the proposed methodology, and the other houses the simulations used in the article.

## Function Folder
There are seven pieces of code that comprise the implementation of the methodology.
- 1_TestStatistics.R: This R code contains all the test statistics considered in the main paper, including the Wald test statistic with pooled and unpooled variance estiamte, the p-value based on Fisher's exact test, the score (Farrington-Manning) test statistics, and the Bayesian posterior probability of the alternative hypothesis with the Jeffreys prior, the uniform prior, and the logit-normal prior.
- 2_RRInitialization.R: This R code is the realization of Algorithm 1 in the main paper, which is used to generate the rejection region for an arbitrary value of the cutoff.
