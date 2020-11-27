# EpiGP
Epistatic Relationship Matrix Based Genomic Prediction of Phenotypes

[![DOI](https://zenodo.org/badge/218787967.svg)](https://zenodo.org/badge/latestdoi/218787967)

Epistatic relationship matrix based genomic prediction of phenotypes does the phenotype prediction based on univariate and bivariate statistical frameworks for two epistatic models in R:

First, Epistatic Random Regression BLUP (ERRBLUP) model as a full epistatic model including all pairwise SNP interactions. ERRBLUP prediction of phenotypes requires ERRBLUP relationship matrix which should be initially calculated.  

Second, selective Epistatic Random Regression BLUP (sERRBLUP) model as a reduced epistatic model. sERRBLUP model contains a desired proportion of pairwise SNP interactions based on the estimated effect sizes or estimated effect variances. sERRBLUP prediction of phenotypes requires sERRBLUP relationship matrix which should be initially calculated based on Pairwise SNP interaction effect sizes or their variances.  

This package provides all the required functions in both univariate and bivariate statistical frameworks step by step, in addition to one function for ERRBLUP and one function for sERRBLUP which do all the analysis in one step.  Besides, the best proportion of pairwise SNP interactions which results in the highest prediction accuracy could be fined by a defined function in the package.

The provided dataset for the example of each function of the package is the wheat dataset markers and respective simulated phenotypes. The simulated phenotypes are also provided in the package.

This repository contains our R-package EpiGP and the highly recomented packages (miraculix / RandomFieldsUtils) which speed up 15 times as fast as the regular matrix multiplications on genotype data in R. 


# License
The EpiGP package is licensed under the GPL-3 License. The GNU General Public License Version 3 (GPL-3) is a free, copyleft license for
software and other kinds of works. More details on GPL-3 License is provided in the following link:
https://www.r-project.org/Licenses/GPL-3


# Update
Version 0.3.0
Defining five new functions for ERRBLUP and sERRBLUP epistatic genomic prediction of phenotypes in bivariate statistical framework  
Simulating two correlated phenotypic traits for bivariate epistatic models


