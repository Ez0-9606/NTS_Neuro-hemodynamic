import numpy as np
import matplotlib.pyplot as plt
import itertools as itr
import random
from PYPACKAGE.parula_color import parula
from sklearn.manifold import Isomap



def random_combination(iterable, r):
    "Random selection from itertools.combinations(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(range(n), r))
    return tuple(pool[i] for i in indices)


rat_idx = 0
from PYPACKAGE.paramFunc import getRatNumber
rat_numb = getRatNumber(rat_idx)
from scipy import io
params = io.loadmat('..\\PARAMETER_DATA\\params.mat')
dt = params['DT'][0][rat_idx]
trial = params['TRIAL'][0][rat_idx]
pre = params['PRE'][0][rat_idx]
post = params['POST'][0][rat_idx]
L = post - pre
trialIdx = params['TrialIDX'][0][rat_idx][0]

# trial = trial-1
# trialIdx = trialIdx[1:]

lag = L*(trialIdx[0] -1)
edges = np.arange(0, trial*L + 10e-8, dt) + lag
t = (edges[1:]+edges[0:-1])*0.5 - lag


""" 1. Load data (Cariovascular) """
from PYPACKAGE import cardiovascular_process as cdp
bp_path = '..\\PROCESSED\\BP_DATA\\bp_data_'+ rat_numb +'_v7.mat'
bp_data = cdp.get_bp(bp_path, cdp.time_edge_2_cardio_idx(edges), filt = False)


""" 2. Get spike firing rate """
from PYPACKAGE import spike_process as sp
spike_path = '..\\PROCESSED\\SPIKE_DATA\\spike_data_' + rat_numb + '_v7.mat'
firing_rate = sp.get_firing_rate(file_path_name = spike_path, bin_t = dt, spike_edges = edges, filt = True, window = round(3/dt))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t, firing_rate[:,0])
ax.plot(t, firing_rate[:,1])


from scipy.signal import detrend
for i in range(firing_rate.shape[1]):

    mx1 = np.zeros((trial,1))

    for j in range(trial):
        sidx = np.argmin(abs(t-j*L))
        eidx = np.argmin(abs(t-(j+1)*L))
        if max(firing_rate[sidx:eidx,i]) == 0:
            firing_rate[sidx:eidx,i] = 0
            mx1[j] = 0
        else:
            tmp1 = firing_rate[sidx:eidx,i]
            mx1[j] = max(tmp1)
            # firing_rate[sidx:eidx,i] = (tmp - min(tmp))/(max(tmp)-min(tmp))
            firing_rate[sidx:eidx,i] = tmp1/mx1[j]
            if j == (trial-1):
                firing_rate[-1,i] = firing_rate[-1,i]/mx1[j]
            del tmp1
    firing_rate[:,i] = firing_rate[:,i]*np.mean(mx1)




""" 3. Isomap """
from sklearn.manifold import Isomap
from ripser import ripser as tda
from scipy.io import savemat
# Calculate Isomap 
print('Isomap finish')
subpopnum = round(firing_rate.shape[1]*0.6)

""" Subpopulation analysis """
repeat = 100
save = np.zeros((repeat,), dtype = 'object')
for i in range(0,repeat):
    subpopidx = random_combination(range(firing_rate.shape[1]), subpopnum)
    subFR = firing_rate[:, random_combination(range(firing_rate.shape[1]), subpopnum)]
    isomap = Isomap(n_components = 6)
    sub_isomap = isomap.fit_transform(subFR)

    # H0 & H1
    results = {'h1': []}
    barcodes = tda(sub_isomap, maxdim=1, coeff=2)['dgms']
    results['h1'] = barcodes[1]
    h1 = results['h1']

    # nearest neighbors
    num = 150
    d = isomap.dist_matrix_
    mx = np.max(d, axis = 0)
    mx = np.max(mx, axis = 0)
    mn = np.min(d, axis = 0)
    mn = np.min(mn, axis = 0)

    th = np.arange(num)/(num-1)
    th = th*(mx-mn) + mn
    avg_neigh_num = np.zeros_like(th)
    for j in range(num):
        neigh_num = np.sum(d < th[j], axis = 1) - 1
        avg_neigh_num[j] = np.mean(neigh_num)
    s = int(np.floor(num*0.05))
    e = int(np.floor(num*0.95))

    # Saving
    save[i] = {
    "subpopIdx": subpopidx,
    "firing_rate":subFR,
    "fr_isomap": sub_isomap,
    "h1": h1,
    "th": np.log10(th[s:e]),
    "avg_neigh_num": np.log10(avg_neigh_num[s:e])
    }
    print(i)
    
savemat("..\\POSTPROCESSED_MAT_CODE\\DATA\\subpopulationAnalysis_matlab_save_60%.mat", {'subpopulationAnalysis':save})


