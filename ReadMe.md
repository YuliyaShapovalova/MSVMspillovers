# MSVMspillovers

Measuring and Quantifying Uncertainty in Volatility Spillovers: A Bayesian Approach

## Citation

This is supplementary code for the 2022 paper:

@article{shapovalova2022spillovers,

  title={Measuring and Quantifying Uncertainty in Volatility Spillovers: A Bayesian Approach},
  
  author={Shapovalova, Yuliya and Eichler, Michael},
  
  journal={arxiv},
  
  year={2021}
}

## Running MCMC

run_MCMC.R -- runs one of the represenative examples from the paper; SMC part depends on the Fortran code -- likfort2020.f90 (see utils); you might have to recompile it to run on you machine to run the code (check Makefile). 

## Computing spillover measures

spillovers_summary.R -- computes spillovers summaries based on the MCMC samples.

