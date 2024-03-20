# Network analysis

After running our MCMC, we can analyse our parameter sets. These are provided in the `param_sets` folder.

We can further analyse this parameter set: `01_MAP_network_analysis.ipynb` provides code to reprouduce panels in Figure 3 of the manuscript.

`02_all_network_analysis.ipynb` is the code to reproduce Figure 4 of the manuscript.

And `03_signalling_perturbations.py` will reproduce the signalling perturbations for Figure 5.

# param_sets

This folder contains the parameter sets presented in this paper.

`MAP_params.csv`: the best fitting overall parameter value.  To obtain the MAP parameter set, I took every parameter set in the last iteration of MCMC, and computed the likelihood score for that parameter set for *all AGETs*. The MAP parameter set is the one with the best score.

`all_networks.csv`: Well-fitting networks filtered as described in the manuscript and presented in commented-out code in `02_all_networks_analysis.ipynb`. The raw output files from MCMC are too large to include in this github repository.

`median_networks_for_clusters.csv`: I clustered the well fitting networks as described in `02_all_network_analysis.ipynb`, then took the median value for each parameter set and used this as representatives of each cluster.

