import numpy as np
import scipy.signal as sig

def bp_hr_preprocess(bp_raw, fs, bp_lfilt_cutoff, hr_window):
    
    bp_preprocess = bp_raw
 
    hr_preprocess = hr_cal(bp_preprocess[:,0], bp_preprocess[:,1], hr_window, fs)

    LF_num, LF_den = sig.butter(2, bp_lfilt_cutoff/(0.5*fs), btype = 'low')
    bp_preprocess[:,1] = sig.filtfilt(LF_num, LF_den, bp_preprocess[:,1])

    return bp_preprocess, hr_preprocess


def hr_cal(bp_t, bp_raw, window, fs, cutoff = 1):
    HF_num, HF_den = sig.butter(2, cutoff/(fs*0.5), 'high')    
    filt = sig.filtfilt(HF_num, HF_den, bp_raw)
    
    dist = np.ceil(0.1*fs)
    peaks, _ = sig.find_peaks(filt, distance = dist)
    
    if peaks[-1] < (len(filt)-1):
        fake_element = len(filt)-1
        peaks = np.append(peaks, fake_element)
        beat_counts = np.bincount(peaks)
        beat_counts[-1] = 0
    else:
        beat_counts = np.bincount(peaks)
    
    # single-beat (instantaneous) heart rate
    hr = beat_counts*fs*60

    # moving average filtering with padding
    padding_front = np.average(hr[0:window])*np.ones(window)
    padding_end = np.average(np.append(hr[-window:-1], hr[-1]))*np.ones(window)    
    moving_avg = 1/window*np.ones(window)
    mv_hr = np.append(padding_front, hr)
    mv_hr = np.append(mv_hr, padding_end)
    mv_hr = sig.filtfilt(moving_avg, 1, mv_hr)
    mv_hr = mv_hr[window:-window]

    #import matplotlib.pyplot as plt
    #plt.figure
    #plt.plot(bp_t, filt)
    #plt.scatter(bp_t[peaks], filt[peaks], c = 'r')
    #plt.plot(mv_hr)
    #plt.show()

    return mv_hr