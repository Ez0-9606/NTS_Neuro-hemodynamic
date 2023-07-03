import numpy as np
import math
import time
import matplotlib.pyplot as plt
from scipy import signal as sig





# Ref) Nat. Neurosci. 19
def gaussian_wind_fn(mu, sigma, x):
    '''Normalized Gaussian'''
    return np.exp(-((x - mu)**2) / (2 * sigma**2)) / (sigma * np.sqrt(2 * np.pi))

def get_kernel_sum(spike_list, t_points, win_fun):
    '''Sum a bunch of kernels centered at each spike time. This can be slow, but 
    is only run once at the beginning so not optimizing.
    '''
    result = np.zeros_like(t_points)
    for spike in spike_list:
        result = result + win_fun(spike, t_points)
    return result





'''
Zero-mean/normalized/zero-phase gaussian filter
'''
def gaussian_filter(L, sigma, bin_t, data):
    x = bin_t*(np.arange(L) - np.floor(L*0.5))
    gaussian_win = np.exp(-(x**2) / (2 * sigma**2)) / (sigma * np.sqrt(2 * np.pi))
    
    return sig.filtfilt(gaussian_win, 1, data)


def get_firing_rate(spike_times, bin_t = 0.5, spike_edges = np.array([])):
    
    n_ch = spike_times.shape
    n_ch = n_ch[0]

    if spike_edges.size == 0:
        global_max = float(0)
        # Find the bin_edges
        for i in range(n_ch):
            local_max = max(spike_times[i])
            if global_max < local_max:
                global_max = local_max
        bin_n = math.ceil(global_max/bin_t)
        #bin_edges = np.arange(0, (bin_n+1)*bin_t, bin_t)
        #time_axis = bin_edges[:-1] + (bin_t / 2.)

    bin_counts = []
    for i in range(n_ch):
        if spike_edges.size == 0:
            temp = spike_times[i] + (bin_n*bin_t,) #add a fake element at the end to make a boundary
            temp = [x/bin_t for x in temp] # (1.2, 4.2, 7.5, ...) -> (0.3, 1.4, 2.5, ...) : range transform from (bin_t)*i <= x < (bin_t)*(i+1) to i <= x <i+1
            temp = np.floor(temp) # (0.3, 1.4, 2.5, ...) -> (0, 1, 2, ...) : make elements the boundary integers for counting
            temp = temp.astype('int') # data type transform from float to the integer for bin-counting
            counts = np.bincount(temp) #(0, 1, 2) => (1 for 0<=x<1, 1 for 1<=x<2, 1 for 2<=x<3, ...) : Counting
            counts = counts[:-1] # remove the fake element
            counts = [x/bin_t for x in counts] #fransform unit into Hz (Samples/second)
            if i == n_ch-1:
                original = counts
        else:
            s_pad = (spike_edges[0]-1,)
            e_pad = (spike_edges[-1]+1,)
            temp = np.digitize(s_pad+spike_times[i]+e_pad, spike_edges, right = False)
            counts = np.bincount(temp)[1:-1]

        L = 11
        sigma = bin_t
        # Apply Gaussian filtering
        counts = gaussian_filter(L, sigma, bin_t, counts)

        bin_counts.append(counts)
        del temp, counts

    return np.array(bin_counts, dtype = 'object')




