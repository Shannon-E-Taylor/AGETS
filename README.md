# Approximated Gene Expression Trajectories (AGETs) for Gene Regulatory Network Inference on Cell Tracks 

This github repository contains code associated with Spiess, Taylor et al. 2022: Approximated Gene Expression Trajectories (AGETs) for Gene Regulatory Network Inference on Cell Tracks.  (TODO add DOI when paper is published). 

Creating Approximated Gene Expression Trajectories (AGETs).ipynb : 
This is the file that created the AGETs we used for the manuscript. 

MCMC with AGETs.ipynb 
This file describes the process for MCMC in our formulation. Note that this code is NOT speed-optimized. 

run_mcmc_fast.py 
This file is the set of scripts that we used for MCMC. It is optimized as far as I could for speed and as such is ~50x faster than the jupyter notebook. 

signalling_perturbations.py 
For reproducibility, this is how Figure 5 of the manuscript was generated. 
