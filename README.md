# MCPF
An R Package Focusing on Association Analysis in the Distance-based Regression Model

## Description
Distance-based regression model has become a powerful approach to identifying  phenotypic associations in many fields. It is found to be particularly useful for high-dimensional biological and genetic data with proper distance or similarity measures being available.
The pseudo F statistic used in this model accumulates information and is effective when the signals, i.e., the variations represented by the eigenvalues of the similarity matrix, scatter evenly along the eigenvectors of the similarity matrix. However, it might lose power for the uneven signals. To deal with this issue, we propose a group analysis on the variations of signals along the eigenvalues of the similarity matrix and take the maximum among them. The new procedure can automatically choose an optimal grouping point on some given thresholds and thus can improve the power evidence. Extensive computer simulations and applications to a prostate cancer data and an aging human brain data illustrate the effectiveness of the proposed method.

The package MCPF contains two functions for the p-value calculation in the distance-based regression model, including simreg (the pseudo F test statistic, abbreviated as PF) and mcpf (maximum of the PF statistics, abbreviated as MCPF). This package also contains a function (calcKerMat_Cen) which can calculate the centered similarity matrix for a given data matrix.

## Installation
The R package MCPF can be installed through the following R code.

> install.packages("devtools")
> library(devtools)
> install_github("sunflowerJJ/MCPF")