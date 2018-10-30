# Nonparametric-MIDAS-with-Fourier-Transformation

## Pre-requists
The code is developed under R (https://www.r-project.org) as well as C++ (See Reference).
Packages which should be installed in the code, have been included in sample.R and functions.R.

## Files Comments
functions.R -- All functions used to do clustering.
sample.R -- An example to use Nonparametric MIDAS with Fourier Transformation
data.csv -- sample data
        15 subjects in each panel; 2 panels in total.
        Sample size = 100.
        Frequency Ratio = 40.
admmmcp.cpp -- ADMM with MCP algorithm.

## Reference
admmmcp.cpp is the ADMM with MCP code in Zhu&Qu 2015. The original code can be found via the following link: https://github.com/Xiaolu-Zhu/LongitudinalClustering.
