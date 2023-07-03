import numpy as np
import math
from scipy.signal import find_peaks
from scipy.signal import butter, filtfilt

#function 'thres_cross' returns the timing array of spikes
#
#(input)
#    filt_raw: 300-5000 Hz high-pass filtered extracellular recording signals [number of channels] x [filtered recording]
#    fs: sampling rate (24414.0625 for Tucker-Davis Technology(TDT) single-unit recording setup)
#    t_bin: bin size of temporal(time) subset to calculate RMS and STD
#    thres: threshold factor which determines the spike detection threshold of RMS + STD*(thres)
#
#(output)
#    spikes: timing data of neural spikes
#
#copyright @ Jiho Lee (2022) from Innovative Medical Solution Laboratory, POSTECH    

def single_unit_filtering(non_filt_raw, fs, LF_cut_freq, HF_cut_freq, order = 3):
    lf_cut = LF_cut_freq / (fs*0.5)
    hf_cut = HF_cut_freq / (fs*0.5)
    bp_num, bp_den = butter(order, [lf_cut, hf_cut], btype = 'band')
    filt_data = filtfilt(bp_num, bp_den, non_filt_raw)
    
    return filt_data



def thres_cross(filt_raw, fs, t_bin, thres):
    
    dt = 1/fs

    logic = filt_raw < 0
    logic = logic.astype(int)

    tst = filt_raw * logic
    
    del logic
    
    data_size = np.shape(filt_raw)
    numb_bin = math.floor(data_size[1]/(fs*t_bin))
    size_bin = math.floor(data_size[1]/numb_bin)

    spikes = []

    for i in range(0, data_size[0]):

        # local rms & std are slightly different with MatLab's results which look due to numerical error
        # however, the overall tendencies for spike trains from Python and Matlab are same
        raw_temp = np.reshape(filt_raw[i, 0:numb_bin*size_bin], (numb_bin, size_bin))
        
        std = np.std(raw_temp, axis = 1)

        raw_temp = raw_temp*raw_temp
        raw_temp = raw_temp.sum(axis = 1)/size_bin
        
        rms = np.sqrt(raw_temp)
        
        del raw_temp

        spike_temp = tuple()
        
        for j in range(0, numb_bin):
            #if j == 1022:
            #    plt.plot(-1*tst[i, j*size_bin:(j+1)*size_bin])
            #    plt.axhline(y=rms[j] + thres*std[j], color='r', linestyle='-')
            #    plt.show()
            #    print('wait')

            peaks, properties = find_peaks(-1*tst[i, j*size_bin:(j+1)*size_bin], height = rms[j] + thres*std[j])
            
            #plt.close()

            #idx = np.append(idx, peaks + j*size_bin)
            spike_temp = spike_temp + tuple((peaks + j*size_bin)*dt)
            del peaks, properties

        spikes.append(spike_temp)
        del spike_temp

    return spikes

def thres_cross_tertrode(filt_raw, fs, t_bin, thres):
    n_tet = filt_raw.shape[0]/4
    spikes = []
    
    for i in range(int(np.ceil(n_tet))):
        if 4*(i+1) > filt_raw.shape[0]:
            temp_tet_spike = thres_cross(filt_raw[4*i:,:], fs, t_bin, thres)
        else:
            temp_tet_spike = thres_cross(filt_raw[4*i:4*(i+1),:], fs, t_bin, thres)

    spikes = spikes + temp_tet_spike
    
    return spikes
