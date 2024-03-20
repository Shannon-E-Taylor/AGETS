import pickle

import matplotlib.pyplot as plt 
from matplotlib import colors

import numpy as np 

import os 
from numba import jit, cfunc, njit
from numbalsoda import lsoda, lsoda_sig

import pandas as pd

from joblib import Parallel, delayed




#######
# scripts to run simulation 
#######

#region 

fastmath = False

@njit(fastmath=fastmath)
def dot(x, y):
    s = 0
    for i in range(len(x)):
        s += x[i]*y[i]
    return s

@njit(fastmath=fastmath)
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
        dot(W3, s) + dot(E3, B) + h[2]
    ]

    # we assign results to du rather than returning them
    du[0] = R[0] * g(u[0]) - lmd[0] * s[0]
    du[1] = R[1] * g(u[1]) - lmd[1] * s[1]
    du[2] = R[2] * g(u[2]) - lmd[2] * s[2]





funcptr = PSH_lsoda.address # address to ODE function

@njit(fastmath=fastmath)
def solve_lsoda(s, t_interval, params): 
    # define jitted functions 
    usol, success = lsoda(funcptr, s, t_interval, params, 
        rtol = 1.0e-8, atol = 1.0e-8) # these are important for acceptable accuracy 
    assert success 
    return usol 



class Simulate_on_tracks_highwnt:
    '''
    The simulator takes as input a set of 24 parameters describing the GRN 
    and uses ODEs with the initial and boundary conditions from the AGETs to 
    simulate gene expression for the chosen 10 AGETs. The simulated gene expression
    is then added to the AGET, so that the AGET includes the original target as well 
    as the simulated gene expression.
    '''

    def __init__(self, parameters, list_of_tracks):
        self.params = parameters
        self.list_of_tracks = list_of_tracks 

    def simulation(self):
        params = self.params 
        tracks_df = self.list_of_tracks

        # precompute t_intervals to speed up function 
        base_t_intervals = []
        for i in range(1, 61): # number of timepoints 
            t_interval = np.linspace(i*90/3600/3, (i+1)*90/3600/3, 10)
            base_t_intervals.append(t_interval)

        # precompute the Time_nPSM
        time_npsm = [i*90/3600/3 for i in range(1, 62)]

        
        # add the params to an empty array of shape (26) to pass to PSH_lsoda 
        params_to_pass = np.empty((26, ))
        params_to_pass[0:24] = params
        
        df_out = [] 
        # Iterate through celltracks (AGETs)
        for df_celltrack in tracks_df:
            # add biological time for the ODEs 
            # I precomputed this earlier to speed up the iterations 
            if df_celltrack.shape[0] == 61: 
                df_celltrack['Time_nPSM'] = time_npsm
            else: 
                df_celltrack['Time_nPSM'] = df_celltrack['Time']*90/3600/3 

            # collect parameters / ICs 
            simulated_expression = [df_celltrack.loc[0, ['g1', 'g2', 'g3']].values]
            B1 = df_celltrack['Wnt'].values
            # B1 = 0.01
            B2 = df_celltrack['FGF'].values
            Time_nPSM = df_celltrack['Time_nPSM'].values
            time = df_celltrack['Time'].values

            # loop through positions in dataframe 
            for index in df_celltrack.index[:-1]: 
                # if the timestep = 1, 
                # we can use the precomputed information above 
                count_timesteps_between = time[index+1] - time[index]
                if count_timesteps_between != 1: 
                    t_interval = np.linspace(Time_nPSM[index], Time_nPSM[index+1], 10*int(count_timesteps_between))
                else: 
                    t_interval = base_t_intervals[int(time[index])-1]

                # define initial conditions 
                # and signalling 
                s0 = simulated_expression[-1]
                params_to_pass[24:26] =  [1.5 ,B2[index]]

                simulated_expression.append(solve_lsoda(
                    s0, 
                    t_interval, 
                    params_to_pass
                    )[-1]
                )

            simulated_expression = np.array(simulated_expression)
            df_celltrack['g1_sim'] = simulated_expression[:, 0]
            df_celltrack['g2_sim'] = simulated_expression[:, 1]
            df_celltrack['g3_sim'] = simulated_expression[:, 2]
            df_out.append(df_celltrack) # append AGET with the simulated expression to the output list

        return df_out




class Simulate_on_tracks_lowfgf:
    '''
    The simulator takes as input a set of 24 parameters describing the GRN 
    and uses ODEs with the initial and boundary conditions from the AGETs to 
    simulate gene expression for the chosen 10 AGETs. The simulated gene expression
    is then added to the AGET, so that the AGET includes the original target as well 
    as the simulated gene expression.
    '''

    def __init__(self, parameters, list_of_tracks):
        self.params = parameters
        self.list_of_tracks = list_of_tracks 

    def simulation(self):
        params = self.params 
        tracks_df = self.list_of_tracks

        # precompute t_intervals to speed up function 
        base_t_intervals = []
        for i in range(1, 61): # number of timepoints 
            t_interval = np.linspace(i*90/3600/3, (i+1)*90/3600/3, 10)
            base_t_intervals.append(t_interval)

        # precompute the Time_nPSM
        time_npsm = [i*90/3600/3 for i in range(1, 62)]

        
        # add the params to an empty array of shape (26) to pass to PSH_lsoda 
        params_to_pass = np.empty((26, ))
        params_to_pass[0:24] = params
        
        df_out = [] 
        # Iterate through celltracks (AGETs)
        for df_celltrack in tracks_df:
            # add biological time for the ODEs 
            # I precomputed this earlier to speed up the iterations 
            if df_celltrack.shape[0] == 61: 
                df_celltrack['Time_nPSM'] = time_npsm
            else: 
                df_celltrack['Time_nPSM'] = df_celltrack['Time']*90/3600/3 

            # collect parameters / ICs 
            simulated_expression = [df_celltrack.loc[0, ['g1', 'g2', 'g3']].values]
            B1 = df_celltrack['Wnt'].values
            # B1 = 0.01
            B2 = df_celltrack['FGF'].values
            Time_nPSM = df_celltrack['Time_nPSM'].values
            time = df_celltrack['Time'].values

            # loop through positions in dataframe 
            for index in df_celltrack.index[:-1]: 
                # if the timestep = 1, 
                # we can use the precomputed information above 
                count_timesteps_between = time[index+1] - time[index]
                if count_timesteps_between != 1: 
                    t_interval = np.linspace(Time_nPSM[index], Time_nPSM[index+1], 10*int(count_timesteps_between))
                else: 
                    t_interval = base_t_intervals[int(time[index])-1]

                # define initial conditions 
                # and signalling 
                s0 = simulated_expression[-1]
                params_to_pass[24:26] =  [B1[index],0.01]

                simulated_expression.append(solve_lsoda(
                    s0, 
                    t_interval, 
                    params_to_pass
                    )[-1]
                )

            simulated_expression = np.array(simulated_expression)
            df_celltrack['g1_sim'] = simulated_expression[:, 0]
            df_celltrack['g2_sim'] = simulated_expression[:, 1]
            df_celltrack['g3_sim'] = simulated_expression[:, 2]
            df_out.append(df_celltrack) # append AGET with the simulated expression to the output list

        return df_out



class Simulate_on_tracks:
    '''
    The simulator takes as input a set of 24 parameters describing the GRN 
    and uses ODEs with the initial and boundary conditions from the AGETs to 
    simulate gene expression for the chosen 10 AGETs. The simulated gene expression
    is then added to the AGET, so that the AGET includes the original target as well 
    as the simulated gene expression.
    '''

    def __init__(self, parameters, list_of_tracks):
        self.params = parameters
        self.list_of_tracks = list_of_tracks 

    def simulation(self):
        params = self.params 
        tracks_df = self.list_of_tracks

        # precompute t_intervals to speed up function 
        base_t_intervals = []
        for i in range(1, 61): # number of timepoints 
            t_interval = np.linspace(i*90/3600/3, (i+1)*90/3600/3, 10)
            base_t_intervals.append(t_interval)

        # precompute the Time_nPSM
        time_npsm = [i*90/3600/3 for i in range(1, 62)]

        
        # add the params to an empty array of shape (26) to pass to PSH_lsoda 
        params_to_pass = np.empty((26, ))
        params_to_pass[0:24] = params
        
        df_out = [] 
        # Iterate through celltracks (AGETs)
        for df_celltrack in tracks_df:
            # add biological time for the ODEs 
            # I precomputed this earlier to speed up the iterations 
            if df_celltrack.shape[0] == 61: 
                df_celltrack['Time_nPSM'] = time_npsm
            else: 
                df_celltrack['Time_nPSM'] = df_celltrack['Time']*90/3600/3 

            # collect parameters / ICs 
            simulated_expression = [df_celltrack.loc[0, ['g1', 'g2', 'g3']].values]
            B1 = df_celltrack['Wnt'].values
            # B1 = 0.01
            B2 = df_celltrack['FGF'].values
            Time_nPSM = df_celltrack['Time_nPSM'].values
            time = df_celltrack['Time'].values

            # loop through positions in dataframe 
            for index in df_celltrack.index[:-1]: 
                # if the timestep = 1, 
                # we can use the precomputed information above 
                count_timesteps_between = time[index+1] - time[index]
                if count_timesteps_between != 1: 
                    t_interval = np.linspace(Time_nPSM[index], Time_nPSM[index+1], 10*int(count_timesteps_between))
                else: 
                    t_interval = base_t_intervals[int(time[index])-1]

                # define initial conditions 
                # and signalling 
                s0 = simulated_expression[-1]
                params_to_pass[24:26] =  [B1[index],B2[index]]

                simulated_expression.append(solve_lsoda(
                    s0, 
                    t_interval, 
                    params_to_pass
                    )[-1]
                )

            simulated_expression = np.array(simulated_expression)
            df_celltrack['g1_sim'] = simulated_expression[:, 0]
            df_celltrack['g2_sim'] = simulated_expression[:, 1]
            df_celltrack['g3_sim'] = simulated_expression[:, 2]
            df_out.append(df_celltrack) # append AGET with the simulated expression to the output list

        return df_out


#endregion 

def run_sim(params): 

    trajectories_tbxta = []
    trajectories_tbx16 = []
    trajectories_tbx24 = []

    #for every cell that made the cut
    for i in range(data_subset.shape[0]):
        # g1, g2, g3 are at these indices 
        IC = np.array(data_subset.loc[i, ['g1', 'g2', 'g3']])

        params_to_pass = np.empty((26, ))
        params_to_pass[0:24] = params
        params_to_pass[24:26] =  data_subset.loc[i, ['Wnt', 'FGF']].values

        single_cell_trajectory = solve_lsoda(IC, t_interval, params_to_pass)
        trajectories_tbxta.append(single_cell_trajectory[:,0])
        trajectories_tbx16.append(single_cell_trajectory[:,1])
        trajectories_tbx24.append(single_cell_trajectory[:,2])
        
    return(
        [np.array(trajectories_tbxta), 
        np.array(trajectories_tbx16), 
        np.array(trajectories_tbx24)])




sim = pd.read_csv(f'../Input/network7_sim_data.csv')


sim_t1 = sim[sim['Time'] == 1]
min_x = min(sim_t1['X'])
max_x = max(sim_t1['X'])

sim_t1['x_norm'] = (sim_t1['X'] - min_x) / (max_x - min_x)
data_subset = sim_t1[sim_t1['x_norm'] < 0.2]
data_subset = data_subset.reset_index(drop = True)
# data_subset['Wnt'] = 0.01

t_interval = np.linspace(0, 200, 10100)

with open("../Input/List_of_all_cell_tracks_starttoend.txt", "rb") as fp:   
    list_of_cell_tracks = pickle.load(fp) # Load the chosen AGETs#

modified_tracks = []

list_of_cell_tracks = [track for track in list_of_cell_tracks if list(track['X'])[-1] < 130]

# sim = pd.concat(Simulate_on_tracks_highwnt(params, list_of_cell_tracks).simulation())




# files = os.listdir('../output/final_t_samples/')
# files = [f for f in files if '40000' in f]

files = np.load('../analysis/good_param_sets.npy', allow_pickle=True)

def run_sims(params, f, num): 
    sim_wnt = pd.concat(Simulate_on_tracks_highwnt(params, list_of_cell_tracks).simulation())
    sim_fgf = pd.concat(Simulate_on_tracks_lowfgf(params, list_of_cell_tracks).simulation())
    sim_wt = pd.concat(Simulate_on_tracks(params, list_of_cell_tracks).simulation())


    sim_wt = sim_wt[sim_wt['Time'] == 61]
    sim_fgf = sim_fgf[sim_fgf['Time'] == 61]
    sim_wnt = sim_wnt[sim_wnt['Time'] == 61]


    for i in range(1, 4): 
        sim_wt[f'g{i}_fgf'] = sim_fgf[f'g{i}_sim']
        sim_wt[f'g{i}_wnt'] = sim_wnt[f'g{i}_sim']

    # return(sim_wt)
    sim_wt.to_csv(f'../output/simulations/simulations_perturbed_signalling{f}_{num}.csv')


files = np.unique(files)
print(files)
print(len(files))



for f in files: 
    # print(j)
    # if True: # overwrite now we're doing the last n parameters 
    if f'../output/simulations/simulations_perturbed_signalling{f}_0.csv' not in os.listdir('../output/simulations/'): 
        all_params = np.loadtxt('../output/final_t_samples/' + f + '.txt_final_timepoint_samp.txt')
        # all_params = all_params[-2]
        Parallel(
            n_jobs=15, backend="multiprocessing"
            )(delayed(run_sims)(params, f, num) for num, params in enumerate(all_params))
        # print(outputs.shape)
        print(all_params.shape)
    else: 
        print ('already done')


