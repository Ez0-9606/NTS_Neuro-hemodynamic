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

def getClosestIdx(target, source):
    "Get idice array corresponding to the values in array 'source' that are closest to each value in 'target'."
    l = len(target)
    idx = np.zeros(l, dtype = int)
    for i in range(0,l):
        idx[i] = np.argmin(abs(source-target[i]))
    return idx

def trialFold(repData, trial):
    "Fold a serial data array of repeated trials into matrix form (Time) x (Repeated Trials)"
    rowL = int(np.floor(repData.shape[0]/float(trial)))
    colL = trial
    foldData = np.zeros((rowL, colL))
    for i in range(0,colL):
        foldData[:,i] = repData[i*rowL:(i+1)*rowL]
    return foldData

def trialFoldMean(repData, trial):
    'repData should be (Time Length x Trial) x (Subject number) matrix'
    tmp = trialFold(repData[:,0],trial)
    foldMean = np.zeros((tmp.shape[0], repData.shape[1]))
    foldMean[:,0] = np.mean(tmp,axis=1)
    for i in range(1, repData.shape[1]):
        tmp = trialFold(repData[:,i], trial)
        foldMean[:,i] = np.mean(tmp,axis=1)
    return foldMean



from PYPACKAGE.paramFunc import getRatNumber
from scipy import io
from PYPACKAGE import cardiovascular_process as cdp
from PYPACKAGE import spike_process as sp

tNorm = np.arange(0, 241)*0.5 - 20
frAll = np.array([])
for rat_idx in range(0,10):
    rat_numb = getRatNumber(rat_idx)
    params = io.loadmat('..\\PARAMETER_DATA\\params.mat')
    dt = params['DT'][0][rat_idx]
    # trial = params['TRIAL'][0][rat_idx]
    trial = 3
    pre = params['PRE'][0][rat_idx]
    post = params['POST'][0][rat_idx]
    L = post - pre

    trialIdx = params['TrialIDX'][0][rat_idx][0]
    if rat_idx == 0:
        trialIdx = trialIdx[1:]

    lag = L*(trialIdx[0] -1)
    edges = np.arange(0, trial*L + 10e-8, dt) + lag
    t = (edges[1:]+edges[0:-1])*0.5 - lag

    """ 1. Load data (Cariovascular) """
    bp_path = '..\\PROCESSED\\BP_DATA\\bp_data_'+ rat_numb +'_v7.mat'
    bp_data = cdp.get_bp(bp_path, cdp.time_edge_2_cardio_idx(edges), filt = False)


    """ 2. Get spike firing rate """
    spike_path = '..\\PROCESSED\\SPIKE_DATA\\spike_data_' + rat_numb + '_v7.mat'
    firing_rate = sp.get_firing_rate(file_path_name = spike_path, bin_t = dt, spike_edges = edges, filt = True, window = round(3/dt))

    from scipy.signal import detrend
    rat_mx = 0
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
        neuron_mx = np.mean(mx1)        
        firing_rate[:,i] = firing_rate[:,i]*neuron_mx
        if rat_mx < neuron_mx:
            rat_mx = neuron_mx
    firing_rate = firing_rate/rat_mx        

    tFold = np.arange(firing_rate.shape[0])*dt + pre
    tNormAll = np.append(tNorm, tNorm+L)
    tNormAll = np.append(tNormAll, tNorm+2*L)

    if rat_idx == 0:
        frAll = firing_rate[getClosestIdx(tNormAll,tFold),:]
        bpAll = bp_data[getClosestIdx(tNormAll,tFold)].reshape(len(tNormAll),1)
    else:
        frAll = np.append(frAll, firing_rate[getClosestIdx(tNormAll,tFold),:], axis=1)
        bpAll = np.append(bpAll, bp_data[getClosestIdx(tNormAll,tFold)].reshape(len(tNormAll),1), axis=1)

    



""" 3. Isomap """
from sklearn.manifold import Isomap
isomap = Isomap(n_components = 6)
# Calculate Isomap 
iso = isomap.fit_transform(frAll)
print('Isomap finish')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(iso[:,0], iso[:,1], c = bp_data[getClosestIdx(tNormAll,tFold)], cmap = 'viridis')




""" 4. Persistent homology """
from ripser import ripser as tda
# H0 & H1
results = {'h1': []}
barcodes = tda(iso, maxdim=1, coeff=2)['dgms']
h1 = barcodes[1]

# H2. Need to subsample points for computational tractability if 
# number of points is large (can go higher but very slow)
if iso.shape[0] > 1500:
    idx = np.random.choice(np.arange(iso.shape[0]), 1500, replace=False)
    H2_rates = iso[idx,:]
else:
    H2_rates = iso
print('h2 check finish')

barcodes = tda(H2_rates, maxdim=2, coeff=2)['dgms']
h2 =  barcodes[2]





""" 5. neighbors vs. deodesic distance """
num = 150
lowDimIsomap = Isomap(n_components = 6)
lowDimIsomap.fit(iso)
d = lowDimIsomap.dist_matrix_
mx = np.max(d, axis = 0)
mx = np.max(mx, axis = 0)
mn = np.min(d, axis = 0)
mn = np.min(mn, axis = 0)

th = np.arange(num)/(num-1)
th = th*(mx-mn) + mn

corr_integral = np.zeros_like(th)

for i in range(num):
    corr_integral[i] = np.sum(np.sum(d < th[i], axis = 1) - 1)
corr_integral = corr_integral/float(num*num)


fig = plt.figure()
ax = fig.add_subplot(111)
e = int(np.floor(num))
ax.scatter(np.log10(th[1:e]), np.log10(corr_integral[1:e]))

print('distance calculation finish')

# Save struture
overallpop = {
    "bp_data":bpAll,
    "firing_rate":frAll,
    "fr_isomap": iso,
    "h1": h1,
    "h2": h2,
    "th": np.log10(th[1:e]),
    "corr_integral": np.log10(corr_integral[1:e])
}


# Set reference data (rat1)
subpopnum = 28
fr_ref = frAll[:,0:subpopnum]
isomap = Isomap(n_components = 6)
iso_ref = isomap.fit_transform(fr_ref)
refsave = {
    "firing_rate":fr_ref,
    "fr_isomap": iso_ref
}

# Exclude reference data
frAll = frAll[:,subpopnum:]
""" Subpopulation analysis """
repeat = 1000
subsamplepop = np.zeros((repeat,), dtype = 'object')
for i in range(0,repeat):
    subpopidx = random_combination(range(frAll.shape[1]), subpopnum)
    subFR = frAll[:, subpopidx]
    isomap = Isomap(n_components = 6)
    sub_isomap = isomap.fit_transform(subFR)

    # Saving
    subsamplepop[i] = {
    "subpopIdx": tuple(map(lambda i: i + subpopnum, subpopidx)),
    "firing_rate":subFR,
    "fr_isomap": sub_isomap,
    }
    print(i)
subsamplesave = {'reference':refsave, 'test':subsamplepop}

from scipy.io import savemat
savemat("..\\POSTPROCESSED_MAT_CODE\\DATA\\overallPopulationAnalysis_matlab_save.mat", {'overallPopulationAnalysis':overallpop, 'subsamplingAnalysis':subsamplesave})


