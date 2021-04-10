"""
estimation of markov state model happens here...
coreMSM: discretises trajectory with core_traj
k-means: Mostly follows tutorial of pyEMMA for k-means
"""

import matplotlib.pyplot as plt
import numpy as np
import pyemma
from msmtools.analysis import eigenvalues as get_eval


def core_traj(data, delta):
    # discretize trajectory on double-well potential using core MSM approach
    # core sets are defined by (<-, -delta], [delta, ->)
    dtraj = np.empty(data.shape[0], dtype=int)  # initialize discrete trajectory
    dtraj[0] = 0 if data[0] < 0 else 1  # take closest core set as initial state
    for i,x in enumerate(data[1:], 1):
        if(x < -delta):
            dtraj[i] = 0
        elif(x >  delta):
            dtraj[i] = 1
        else:
            dtraj[i] = dtraj[i-1]
    return dtraj

# load trajectories
container = np.load('sim_data/dw_slt_10k_001/20trajs_10000_001.npz')
data = [container[key].reshape(container[key].size, -1) for key in container]

# Settings of MSM estimation
delta = .6
its = []
its2 = []
dt = 0.01
lag_time = 10  
tau = lag_time/dt # pyemma.msm.estimate_markov_model takes tau in simulation timesteps
traj_lens = np.linspace(1000,10000,num=10, dtype=int)


# loop trough single long trajectories
for traj in data:
    #dtraj = core_traj(traj, delta)  # discretise trajectories for coreset approach 
    for t_end in traj_lens:
        # dtraj_selection = dtraj[0:int(t_end/dt)+1]  # selection of core trajectory

        traj_selection  = traj[0:t_end*100+1]
        dtraj_selection = pyemma.coordinates.cluster_kmeans(traj_selection, k=20).dtrajs
        
        # estimate MSM
        try:
            msm = pyemma.msm.estimate_markov_model(dtraj_selection, lag=tau, reversible=True, count_mode='sample')
            
            # store time scale of the model
            msm_timescale_estimate = msm.timescales(1)*dt
            its.append(msm_timescale_estimate if msm_timescale_estimate.size != 0 else np.array([-1]))
        except:
            if np.sum(traj_selection)==0 or np.sum(traj_selection)==int(t_end/dt)+1:
                print("no transition in trajectory, setting timescale to -1")
            else:
                print("something went wrong, setting timescale to -1")
            its.append(np.array([-1]))

