# Constructing AGETs from data

Here is the data and code to construct AGETs from HCR data.
The input data is an excel file exported from Imaris, containing quantifications of gene expression and signalling profiles in the nuclei of different genes.
We also input a cell tracking dataset, formatted as an `.csv`. Please see Speiss and Taylor 2022, Thompson et al. 2021 for details on how these data were generated from microscopy images.

## Folder structre

For details on how to construct the AGETs, please see `Creating Approximated Gene Expression Trajectories (AGETs).ipynb`.
You will need to run this in the provided conda environment.
This code is quite flexible: there are a number of options for troubleshooting provided that we do not discuss in the main manuscript.

`compare_agets.ipynb` will reproduce the graphs in Figure 2 of the manuscript.

### Source images

This folder and subfolders contain the exported Imaris files quantifying gene expression and signalling profiles in nuclei.
Code to access the Imaris files are in the jupyter notebook.

### Tracking data

This folder contains the cell tracks used in this analysis.
