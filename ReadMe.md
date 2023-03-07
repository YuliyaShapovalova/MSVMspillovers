# MSVMspillovers

Measuring and Quantifying Uncertainty in Volatility Spillovers: A Bayesian Approach

## Citation

This is supplementary code for the 2022 paper:

@article{doi:10.1080/26941899.2023.2176379,
author = {Yuliya Shapovalova and Michael Eichler},
title = {Measuring and Quantifying Uncertainty in Volatility Spillovers: A Bayesian Approach},
journal = {Data Science in Science},
volume = {2},
number = {1},
pages = {2176379},
year  = {2023},
publisher = {Taylor & Francis},
doi = {10.1080/26941899.2023.2176379},
URL = {https://doi.org/10.1080/26941899.2023.2176379},
eprint = {https://doi.org/10.1080/26941899.2023.2176379}
}

## Running MCMC

run_MCMC.R -- runs one of the represenative examples from the paper; SMC part depends on the Fortran code -- likfort2020.f90 (see utils); you might have to recompile it to run on you machine to run the code (check Makefile). 

## Computing spillover measures

spillovers_summary.R -- computes spillovers summaries based on the MCMC samples.

