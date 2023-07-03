from scipy import io
import numpy as np
import math
import time
import matplotlib.pyplot as plt
from scipy import signal as sig

def get_spike_timing(file_path_name):
    spike_data = io.loadmat(file_path_name)
    spike_data = spike_data['spike'][:,0]
    spike_times = ()
    for i in range(len(spike_data)):
        add_time = spike_data[i]
        if add_time.shape[0] == 0:
            add_time = ((),)
        else:
            add_time = (tuple(add_time[:,0]),)
        spike_times = spike_times + add_time
    spike_times = np.array(spike_times, dtype = 'object')
    return spike_times

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


#Zero-mean/normalized/zero-phase gaussian filter
def gaussian_filter(L, sigma, bin_t, data):
    x = bin_t*(np.arange(L) - np.floor(L*0.5))
    gaussian_win = np.exp(-(x**2) / (2 * sigma**2)) / (sigma * np.sqrt(2 * np.pi))
    
    return sig.filtfilt(gaussian_win, np.sum(gaussian_win), data)

# remove numerical error of gaussian filtering
def correct_post_filtering_error(firing_rate):
    
    #zero_idx = [(j,i) for i in range(firing_rate.shape[1])
    #    for j, val in enumerate(firing_rate[:,i]) if val < 0]
    zero_idx = [i for i, val in enumerate(np.squeeze(firing_rate)) if val < 0]
    for i in zero_idx:
        # print(i)
        # print(firing_rate[i])
        firing_rate[i] = 0
    return firing_rate

# Attenuate the edge of firing rate 
# From 'Mike Lawrence' answering to a thread in StackExchange/Signal Processing
# https://dsp.stackexchange.com/questions/69643/advice-on-designing-a-digital-filter-that-doesnt-have-phase-sensitive-edge-arti
def attenuate_edges(signal,edge_attenuation_percent):
    time = np.arange(0, len(signal))
    start = int(np.floor(len(time)*edge_attenuation_percent))
    end = int(len(time)-start)
    ramp = (1-np.cos(np.pi*(np.arange(start)/start)))/2
    edge_attenuator = np.ones(len(time))
    edge_attenuator[0:start] = ramp
    edge_attenuator[end:len(time)] = np.flip(ramp)
    return(signal*edge_attenuator)

def edges_padding(signal,edge_attenuation_percent):
    time = np.arange(0, len(signal))
    start = int(np.floor(len(time)*edge_attenuation_percent))
    end = int(len(time)-start)
    left_padd = signal[start:2*start]
    right_padd = signal[end-start:end]
    tmp = signal
    tmp[0:start] = left_padd
    tmp[end:len(time)] = right_padd
    return(tmp)

# normalize firing rate
def firing_rate_normalize(firing_rate):

    for i in range(firing_rate.shape[1]):
        if max(firing_rate[:,i]) < 0.0001:
            continue
        firing_rate[:,i] = firing_rate[:,i]/max(firing_rate[:,i])
    return firing_rate

def get_firing_rate(file_path_name, bin_t = 0.5, spike_edges = np.array([]), normalized = False, filt = False, window = 11):
    
    spike_times = get_spike_timing(file_path_name)
    
    n_ch = spike_times.shape[0]

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

        if filt == True:
            L = window
            sigma = bin_t
            # Apply Gaussian filtering
            counts = gaussian_filter(L, sigma, bin_t, counts)
            counts = correct_post_filtering_error(counts)
            counts = edges_padding(counts, 0.005)
        
        bin_counts.append(counts)
        del temp, counts

    output = np.array(bin_counts, dtype = 'object')
    output = output.transpose()
    output = output.astype('float')

    if normalized == 1:
        output = firing_rate_normalize(output)


    return output
