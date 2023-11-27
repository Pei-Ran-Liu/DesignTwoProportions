# Design of Two Sample Proportions

This repository contains R code associated with the article titled "Are the Tests Overpowered or Underpowered? A Unified Solution to Correctly Specify Type I Errors in the Design of Clinical Trials for Two Sample Proportions." This article was submitted to *Statistics in Medicine*. Please fell free to email peiran.liu@uconn.edu if you have any questions or comments.

The codes aid in the design of two-sample proportion tests in their difference using exact calculations. The repository includes two folders: one contains the functions for the proposed methodology, and the other houses the simulations used in the article.

## Functions Folder
There are seven pieces of code that comprise the implementation of the methodology.
- 1_TestStatistics.R: This R code contains all **the test statistics** considered in the main paper, including the Wald test statistic with pooled and unpooled variance estiamte, the p-value based on Fisher's exact test, the score (Farrington-Manning) test statistics, and the Bayesian posterior probability of the alternative hypothesis with the Jeffreys prior, the uniform prior, and the logit-normal prior.
- 2_RRInitialization.R: This R code is the realization of **Algorithm 1** in the main paper, which is used to generate the boundary of the rejection region for an arbitrary value of the cutoff. The functions are written separately for the frequentist approaches and the Bayesian approaches.
- 3_OCEvaluation.R: This R code contains several functions for calculating **the operating characteristics** of the deterministic decision rule, including the actual type I error, the power, and the value of a test statistic at the boundary of the rejection region.
- 4_TuningProcedure.R: This R code contains four functions corresponding to the realization of **the forward procedure and the backward procedure** proposed in the article for the frequentist approaches and the Bayesian approaches, respectively.
- 5_Deterministic.R: This R code consolidates all the functions mentioned above into a single function. Please note: it's essential to run these five functions in sequence.
- 6_Probabilistic.R: This R code comprises several functions related to the probabilistic component of **the probabilistic test**, or equivalently, the randomized decision rule. It corresponds to Algorithm 3 in the main paper. Please refer to the simulation codes to understand the utilization of these functions.
- 7_GridSearch.R: This R code facilitates **the grid search algorithm** for comparison purposes in measuring run time, showcasing the effectiveness of our proposed method. Please be aware that the execution of this approach might take a considerable amount of time.

## Simulations Folder
The names of these three pieces of code correspond to the superiority trial, the noninferiority trial, and the comparison in running time, as discussed in Section 5 of the main paper. Please refer to the article for detailed settings of these simulation studies.
