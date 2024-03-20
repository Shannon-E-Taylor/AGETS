# Approximated Gene Expression Trajectories (AGETs) for Gene Regulatory Network Inference on Cell Tracks

This github repository contains code associated with Spiess, Taylor et al. 2022: Approximated Gene Expression Trajectories (AGETs) for Gene Regulatory Network Inference on Cell Tracks.  (TODO add DOI when paper is published).

This work is divided up into three sections: 

Step 1: construct the AGETs. The procedure and code for this is described in the folder `01_AGET_construction`. 

Step 2: reverse engineer the GRNs using the constructed AGETs. The code and dependencies for this are in `02_Reverse_engineering_GRNs_from_AGETs`. This will need to be run on a supercomputer. 

Step 3: finally, we validate and investigate our parameter sets. For our system, code for this is presented in  `03_Investigations_using_MAP_network`.

You will need the AGET conda environment

`conda create -n AGET --python=3.6 --file environment.yml`

Then activate with

`conda activate AGET`
