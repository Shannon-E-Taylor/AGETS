#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numpy.random import default_rng
import pandas as pd
import pickle
from scipy.integrate import odeint
from multiprocessing import Pool
from multiprocessing import cpu_count
import time
from numba import jit, cfunc, njit
from numbalsoda import lsoda, lsoda_sig
# for now, this is only installed in the AGETS_ptemcee environment
import emcee
import math
import sys

# info for jit 
fastmath = TRUE


##########
# INPUTS #
##########

# read the number of walkers to run with emcee. 
# increasing the number of walkers can speed up convergence but will slow down the computation 
nwalkers = int(sys.argv[1]) 

# read the number of timesteps (samples) to run emcee for 
# In our system, we needed 10000 samples, but this can vary 
Nsamples = int(sys.argv[2])

# seed for AGET choice  
random_num = int(sys.argv[3])

# the number of AGETs to use for reverse engineering 
# recommend about 10% for start 
n_agets = int(sys.argv[4])

# Increment this so repeated runs don't overwrite previous ones!! 
techrep = int(sys.argv[5])

print('Parameters are: ')
print('nwalkers: ' + str(nwalkers))
print('Nsamples: ' + str(Nsamples))

# name the run based on the paramters above. 
i_run = f'{n_agets}_AGETs.{random_num}_randomnumber.{Nsamples}_it.{nwalkers}_walkers.{techrep}_technicalreplicate'

print('i_run: ' + str(i_run))

###
# Choose AGETs to fit to 
with open("Input/List_of_all_cell_tracks_starttoend.txt", "rb") as fp:
    list_of_cell_tracks = pickle.load(fp) # Load the chosen AGETs#

# we remove anterior tracks 
keep_tracks = [track for track in list_of_cell_tracks if track['X'].values[-1] < 170]

# initialize the random seed, and choose data 
rng = default_rng(random_num)
list_of_cell_tracks=rng.choice(keep_tracks, n_agets, replace = False)



#############
# FUNCTIONS #
#############

# these are the functions to run PSH etc
# they are compiled with jit so that they are fast

# code to run the matrix dot product very quickly 
@njit()
def dot(x, y):
    s = 0
    for i in range(len(x)):
        s += x[i]*y[i]
    return s


@njit()
def g(x):
    return 0.5 * ((x / np.sqrt(x ** 2 + 1)) + 1)

@cfunc(lsoda_sig)
def PSH_lsoda(t, s, du, params):
    # unpack params
    W1 = [params[0], params[1], params[2]]
    W2 = [params[3], params[4], params[5]]
    W3 = [params[6], params[7], params[8]]

    E1 = [params[18], params[17]]
    E2 = [params[19], params[15]]
    E3 = [params[20], params[16]]

    R = [params[9], params[10], params[11]]
    lmd = [params[12], params[13], params[14]]
    h = [params[21], params[22], params[23]]

    B = [params[24], params[25]]

    # do computation
    u = [
        dot(W1, s) + dot(E1, B) + h[0],
        dot(W2, s) + dot(E2, B) + h[1],
        dot(W3, s) + dot(E3, B) + h[2],
    ]

    d_tbxta_dt = R[0] * g(u[0]) - lmd[0] * s[0]
    d_tbx16_dt = R[1] * g(u[1]) - lmd[1] * s[1]
    d_tbx24_dt = R[2] * g(u[2]) - lmd[2] * s[2]

    # we assign results to du rather than returning them
    du[0] = d_tbxta_dt
    du[1] = d_tbx16_dt
    du[2] = d_tbx24_dt


funcptr = PSH_lsoda.address # address to ODE function

@njit()
def solve_lsoda(s, t_interval, params):
    # define jitted functions
    usol, success = lsoda(funcptr, s, t_interval, params,
        rtol = 1.0e-8, atol = 1.0e-8) # these are important for acceptable accuracy
    if success:
        return usol
    else:
        return usol * 0 

print('defined solver')

def precompute_df(tracks_df):
    '''
    This function takes the AGETs identified above, and pre-processes them such that we can run 
    the ODEs with different parameter values quickly on them. 
    '''
    precomputed_information = []

    for df_celltrack in tracks_df:
        df_celltrack = df_celltrack.reset_index(drop = True)

        # get ICs
        s0 = df_celltrack.loc[0, ['g1', 'g2', 'g3']].values

        # extract the gene expression values from the dataframe 
        gene_expression = df_celltrack.loc[:, ['g1', 'g2', 'g3']].values

        # extract the signalling values 
        B1 = df_celltrack['Wnt'].values
        B2 = df_celltrack['FGF'].values

        # Convert microscopy frame number to biological time
        df_celltrack.Time_nPSM = df_celltrack.Time*90/3600/3 

        # loop through dataframe to get time intervals for simulation
        # we simulate the ODE on 10 timepoints between each microscopy frame 
        # each microscopy frame differs by ~0.05AU so we have sufficient temporal resolution 
        t_eval = []
        for index in df_celltrack.index[:-1]:
            t_eval.append(
                np.array(np.linspace(
                    df_celltrack.Time_nPSM[index], df_celltrack.Time_nPSM[index+1], 10)
            ))

        # package the information together 
        precomputed_information.append(
            np.array(
                [np.array(B1),
                np.array(B2),
                np.array(s0),
                np.array(t_eval),
                gene_expression
                ],dtype='object'
                )
        )

    return precomputed_information


@jit(nopython = True, fastmath=fastmath)
def simulate_single_track(B1, B2, t_eval, s0, params):
    '''
    Runs the simulation for the signalling, initial conditions, timesteps, and params
    for one cell track

    B1: Wnt
    B2: FGF
    t_eval: list of timesteps to run the simulation for
    s0: initial conditions
    params: parameter values
    '''

    # define a data object to hold our simulated data #
    # it's faster to do this now, than append to a numpy array later 
    simulated_expression = np.empty((B1.shape[0], 3), dtype='float64')
    simulated_expression[0] = s0

    # iterate through the signalling profile
    for index in range(B1.shape[0]-1):
        s0 = simulated_expression[index]
        params_to_pass = np.concatenate(
            (params, np.array([B1[index], B2[index]])), )
        # run computation
        s1 = solve_lsoda(
            s0,
            t_eval[index],
            params_to_pass
            )[-1]
        simulated_expression[index+1] = s1
    return(simulated_expression)

# @jit(nopython = False, forceobj=True)
def calculate_logl(params, list_of_tracks):
    '''
    Runs the simulation for all tracks in list_of_tracks and params
    And returns the loglikleyhood

    Note that list_of_tracks must be the output from precompute_df()
    '''
    all_simulated_expression = []
    real_expression = []
    for track in list_of_tracks:
        B1, B2, s0, t_eval, geneexp = track
        simulated_expression = simulate_single_track(B1, B2, t_eval, s0, params)
        all_simulated_expression.append(simulated_expression)
        real_expression.append(geneexp)

    # now that we have run the simulation on all cell tracks, 
    # we can concatenate the data together 
    all_simulated_expression = np.concatenate(all_simulated_expression)
    real_expression = np.concatenate(real_expression)

    # now we compute the likelihood 
    # note the sd value differs for each gene 
    ll1 = -0.5 * np.sum(((all_simulated_expression[:, 0] - real_expression[:, 0])/0.2)**2)
    ll2 = -0.5 * np.sum(((all_simulated_expression[:, 1] - real_expression[:, 1])/0.2)**2) 
    ll3 = -0.5 * np.sum(((all_simulated_expression[:, 2] - real_expression[:, 2])/0.1)**2) 
    return ll1 + ll2 + ll3


# wrapper function 
def loglikelihood(theta):
    ll = calculate_logl(theta, precomputed_data)
    if not ll or math.isnan(ll):
        return 0
    else:
        return ll


@njit()
def gaussian(x, mu=0.0, sigma=1.0):
    x = float(x - mu) / sigma
    return math.exp(-x*x/2.0) / math.sqrt(2.0*math.pi) / sigma

def logposterior(theta):
    lp = logprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + calculate_logl(theta, precomputed_data)


precomputed_data = precompute_df(list_of_cell_tracks)



# define the maximum / sigma values for the priors
prior_max = [200 for i in range(8+1+3)] + \
            [200 for i in range(3)] + \
            [200 for i in range(9)]

# these are our positive priors
positive_prior_idx = [
    0, 4, 8, 9, 10, 11, 12, 13, 14, 18
]

# now we ensure that all positive priors are indeed positive 
prior_min = [i * -1 for i in prior_max]
for val in positive_prior_idx:
    prior_min[val] = 0


def logprior(theta): # prior was chosen based on biological intuition for the parameters
    #  check if any of our priors are outside the allowed values
    # if they are return negative infinity 
    # which is a hack to make inappropriate priors return ridiculously poor scores 
    if all([theta[i] < prior_max[i] and theta[i] > prior_min[i] for i in range(24)]):
        lp = 0
    else:
        lp = -np.inf
    return lp



# now we set initial conditions
# with a normal distribution around zero 
# we take the absolute value for positive priors 
p0 = []

for i in range(24):
    init = np.random.normal(0, 1, size=(nwalkers))
    if i in positive_prior_idx:
        init = np.array([abs(val) for val in init])
    p0.append(init)

p0 = np.array(p0).T

ndim = 24 # needed below
 # if running multiple times, i_run is used for the name of the saved results
ncpu = cpu_count() # finde number of available CPUs for parallel processing
print("{0} CPUs".format(ncpu))


# Run MCMC (this takes time)
with Pool(ncpu) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, logposterior, pool=pool)
    start = time.time()
    sampler.run_mcmc(p0, Nsamples, progress=True)
    end = time.time()
    multi_time = end - start
    print("Multiprocessing took {0:.1f} seconds".format(multi_time))


# Save samples
samples = sampler.get_chain(flat=True)
print('samples')
print(samples.shape)
print(type(samples))
print(' ')
with open(f"output/samples_run{i_run}.txt", "wb") as fp:   #Pickling
    pickle.dump(samples, fp)

# Save log probabilities
log_probs = sampler.flatlnprobability
print('log_probs')
print(log_probs.shape)
print(type(log_probs))
print(' ')
with open(f"output/log_probs_run{i_run}.txt", "wb") as fp:   #Pickling
    pickle.dump(log_probs, fp)

# Save MAP params (params with highest posterior probability)
map_params = samples[np.argmax(sampler.flatlnprobability)]
print('map_params')
print(map_params.shape)
print(type(map_params))
print(' ')
with open(f"output/map_params_run{i_run}.txt", "wb") as fp:   #Pickling
    pickle.dump(map_params, fp)

# Save acceptance fraction: MCMC diagnostic
acceptance_fraction = sampler.acceptance_fraction
print('acceptance_fraction')
print(acceptance_fraction.shape)
print(type(acceptance_fraction))
print(' ')
with open(f"output/acceptance_fraction_run{i_run}.txt", "wb") as fp:   #Pickling
    pickle.dump(acceptance_fraction, fp)

# Save autocorrelation time: MCMC diagnostic
autocorr_time = sampler.get_autocorr_time(quiet=True)
print('autocorr_time')
print(autocorr_time.shape)
print(type(autocorr_time))
print(' ')
with open(f"output/autocorr_time_run{i_run}.txt", "wb") as fp:   #Pickling
    pickle.dump(autocorr_time, fp)
