# EpiGP
Epistatic Relationship Matrix Based Genomic Prediction of Phenotypes

EpiGP package does the phenotype prediction based on two statistical models:
First, Epistatic Random Regression BLUP (ERRBLUP) model as a full epistatic model with all pairwise SNP interactions. For ERRBLUP prediction of phenotypes, ERRBLUP relationship matrix should be initially calculated. 
Second, selective Epistatic Random Regression BLUP (sERRBLUP) model as a reduced epistatic model. A selected subset of pairwise SNP interactions based on estimated effects or estimated effects variances are included in sERRBLUP model. For sERRBLUP prediction of phenotypess, sERRBLUP relationship matrix should be initially calculated based on Pairwise SNP interaction effects or their variances. These could be done by the EpiGp package.

This repository contains our R-package EpiGP and the highly recomented packages (miraculix / RandomFieldsUtils) which speed up 15 times as fast as the regular matrix multiplications on genotype data in R. 






[![DOI](https://zenodo.org/badge/218787967.svg)](https://zenodo.org/badge/latestdoi/218787967)
