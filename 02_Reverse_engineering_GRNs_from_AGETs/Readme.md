# Reverse engineering GRNs from AGETs

This folder contains the AGETs and script required to reverse-engineer GRNs.

We have two files that will both run the MCMC:

`MCMC with AGETs.ipynb` describes how the code works and provies a minimal example. This code is more readable.

`run_mcmc_fast.py` has been heavily optimized for speed and uses the numba and numbalsoda libraries to solve the ODEs. I do not believe that further speed optimization is possible in python: if you are wanting to try this method it would be worth considering implementing the code in julia or a similar faster language. 

The code, overall, will: 

1) pre-process and package the AGETs generated in part (1) so that they can run efficiently.

2) Define and precompile the ODEs required for simulation for efficient computation.  

3) Defines the parameters required for MCMC- specifically the priors, initial guesses, and likelhood function.

4) Finally, MCMC is run using the emcee package and the data is saved.  

This code should be run from the terminal with MCMC parameters as arguments, in the following order:

1) number of walkers (eg 100). The `emcee` docs recommend no less than twice the number of cores on your machine, for maximum efficiency, and then as many as you can bear. Experimenting with this number will impact the efficiency of the MCMC but in my experience does not change the sorts of parameters you will get out. More walkers is better but will ~linearly increase the iteration time for MCMC.

2) Number of samples / number of iterations to run MCMC for. We use 10,000 - 20,000 to ensure statistically robustness.

3) Random_num: seed for the random number generator when selecting the AGETs to fit to: useful for reproducibility.

4) Number of AGETs to fit to. This number will require experimentation to ensure your AGETs fully represent your dataset. I  recommend 10% as a starting point but this will be highly dependent on the complexity of your data.

5) technical replicate: Number of times you have run MCMC for this random seed- this is only included so repeated runs do not overwrite your parameter sets!

I recommend running a minimum of five 'biological replicates' (different seeds for the AGET choice, representing different cohorts of AGETs used for reverse engineering), and three 'technical replicates' (repeated runs of MCMC on the same AGET set). This is important to ensure you obtain a range of parameter sets for further analysis.

Example code to run:

`python run_mcmc_fast.py 100 20000 42 200 1`

This code was run in the University of Oxford ARC servers, with one core and 48 (max number of) nodes. Each run took ~12 hours.

# Quick tips on optimization

I spent a lot of time trying to optimize the number of walkers, number of iterations, etc to optimize MCMC efficiency, but I found that by far the most effective change was increasing the number of AGETs I used for reverse engineering. I *strongly recommend* spending a good amount of time validating whether your AGETs adequately represent your biological system of interest before reverse-engineering parameter sets with MCMC. If I was starting this project again, I would also do my inference using Julia, as the language is a lot faster, and there appears to be more flexibility in choosing different MCMC methodologies than is available with emcee. I would also try to reduce the number of parameters as far as possible: any more than 24 would be very difficult with our system.
