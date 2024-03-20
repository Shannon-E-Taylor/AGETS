# Approximated Gene Expression Trajectories (AGETs) for Gene Regulatory Network Inference on Cell Tracks

This github repository contains code associated with Spiess, Taylor et al. 2022: Approximated Gene Expression Trajectories (AGETs) for Gene Regulatory Network Inference on Cell Tracks.  (TODO add DOI when paper is published).

The documents are split into three chunks:

1) Construct the AGETs - see `01_AGET_construction`

2) Reverse engineer the GRNs using these agets: see `02_Reverse_engineering_GRNs_from_AGETs`.

3) Finally, we want to validate and investigate our parameter sets: code for this is in `03_Investigations_using_MAP_network`.

You will need the AGET conda environment

`conda create -n AGET --python=3.6 --file environment.yml`

Then activate with

`conda activate AGET`