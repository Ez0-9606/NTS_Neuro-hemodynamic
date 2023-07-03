import numpy as np
import matplotlib.pyplot as plt
import itertools as itr
import random
from PYPACKAGE.parula_color import parula
from sklearn.manifold import Isomap


# !!!! Only rat 'B15\a1' is calculated for tutorial
#  change the range(X,10) where X = 0~9
for repeat in range(9,10):
    rat_idx = repeat
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

    lag = L*(trialIdx[0] -1)
    edges = np.arange(0, trial*L + 10e-8, dt) + lag
    t = (edges[1:]+edges[0:-1])*0.5 - lag


    """ 1. Load data (Cariovascular) """
    from PYPACKAGE import cardiovascular_process as cdp
    bp_path = '..\\PROCESSED\\BP_DATA\\bp_data_'+ rat_numb +'_v7.mat'
    bp_data = cdp.get_bp(bp_path, cdp.time_edge_2_cardio_idx(edges), filt = False)

    # # plot cardiovascular data
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.scatter(t, bp_data, c = t, cmap = parula())
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.scatter(t, bp_data, c = bp_data, cmap = parula())



    """ 2. Get spike firing rate """
    from PYPACKAGE import spike_process as sp
    spike_path = '..\\PROCESSED\\SPIKE_DATA\\spike_data_' + rat_numb + '_v7.mat'
    firing_rate = sp.get_firing_rate(file_path_name = spike_path, bin_t = dt, spike_edges = edges, filt = True, window = round(3/dt))

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(t, firing_rate[:,0])
    # ax.plot(t, firing_rate[:,1])


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

    mx2 = np.zeros((trial,1))
    for j in range(trial):
        sidx = np.argmin(abs(t-j*L))
        eidx = np.argmin(abs(t-(j+1)*L))
        if min(bp_data[sidx:eidx]) == 0:
                    bp_data[sidx:eidx] = 0
                    mx2[j] = 0
        else:
            tmp2 = bp_data[sidx:eidx]
            mx2[j] = min(bp_data[sidx:sidx+100])
            # firing_rate[sidx:eidx,i] = (tmp - min(tmp))/(max(tmp)-min(tmp))
            bp_data[sidx:eidx] = tmp2/mx2[j]
            if j == (trial-1):
                bp_data[-1] = bp_data[-1]/mx2[j]
            del tmp2
    bp_data = bp_data*np.mean(mx2)

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(t, bp_data)

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(t, firing_rate[:,0])
    # ax.plot(t, firing_rate[:,1])

    print('Data load finish')





    """ 3. Isomap """
    from sklearn.decomposition import PCA
    pca = PCA(n_components = 6, svd_solver='arpack')
    fr_pca = pca.fit_transform(firing_rate)


    from sklearn.manifold import Isomap

    # Calculate Isomap 
    isomap = Isomap(n_components = 6)
    fr_isomap = isomap.fit_transform(firing_rate)
    print('Isomap finish')

    # Isomap visualization
    fig = plt.figure()
    ax_iso_3d = fig.add_subplot(111, projection = '3d')
    ax_iso_3d.scatter(fr_isomap[:,0],fr_isomap[:,1],fr_isomap[:,2], c=bp_data, cmap = parula())

    fig = plt.figure()
    ax_iso_3d = fig.add_subplot(111, projection = '3d')
    ax_iso_3d.scatter(fr_isomap[:,0],fr_isomap[:,1],fr_isomap[:,2], c=t%L, cmap = parula())

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(t, fr_isomap[:,0])
    ax.plot(t, fr_isomap[:,1])





    """ 4. Persistent cohomology """
    from ripser import ripser as tda

    # H0 & H1
    results = {'h0': [], 'h1': [], 'h2': []}
    H1_rates = fr_isomap
    barcodes = tda(H1_rates, maxdim=1, coeff=2)['dgms']
    results['h0'] = barcodes[0]
    results['h1'] = barcodes[1]

    print('h0 h1 finish')


    # H2. Need to subsample points for computational tractability if 
    # number of points is large (can go higher but very slow)
    if fr_isomap.shape[0] > 1500:
        idx = np.random.choice(np.arange(fr_isomap.shape[0]), 1500, replace=False)
        H2_rates = fr_isomap[idx,:]
    else:
        H2_rates = fr_isomap

    print('h2 check finish')

    barcodes = tda(H2_rates, maxdim=2, coeff=2)['dgms']
    results['h2'] = barcodes[2]

    print('h2 finish')


    # Plot Betti barcode H0, H1, H2
    import matplotlib.gridspec as gridspec
    color_list = ['r', 'g', 'm', 'c']
    h0, h1, h2 = results['h0'], results['h1'], results['h2']
    # replace the infinity bar (-1) in H0 by a really large number

    h0[~np.isfinite(h0)] = 100
    # Plot the longest barcodes only
    #plot_prcnt = [99, 98] # order is h0, h1, h2

    # fig = plt.figure()
    # if len(h2) == 0:
    #     plot_prcnt = [50, 50]
    #     to_plot = []
    #     for curr_h, cutoff in zip([h0, h1], plot_prcnt):
    #             bar_lens = curr_h[:,1] - curr_h[:,0]
    #             plot_h = curr_h[bar_lens > np.percentile(bar_lens, cutoff)]
    #             to_plot.append(plot_h)
        
    #     gs = gridspec.GridSpec(2, 1)
    # else:
    #     # plot_prcnt = [98, 92, 80]
    #     plot_prcnt = [25, 25, 25]
    #     to_plot = []
    #     for curr_h, cutoff in zip([h0, h1, h2], plot_prcnt):
    #             bar_lens = curr_h[:,1] - curr_h[:,0]
    #             plot_h = curr_h[bar_lens > np.percentile(bar_lens, cutoff)]
    #             to_plot.append(plot_h)
    #     gs = gridspec.GridSpec(3, 1)

    # for curr_betti, curr_bar in enumerate(to_plot):
    #     ax = fig.add_subplot(gs[curr_betti, 0])
    #     for i, interval in enumerate(reversed(curr_bar)):
    #         ax.plot([interval[0], interval[1]], [i, i], color=color_list[curr_betti], lw=1.5)
    #     # ax.set_xlim([-0.05, 3.05])





    """ 5. neighbors vs. deodesic distance """
    num = 150
    # d = isomap.dist_matrix_
    lowDimIsomap = Isomap(n_components = 6)
    lowDimIsomap.fit(fr_isomap)
    d = lowDimIsomap.dist_matrix_
    # from sklearn.metrics.pairwise import euclidean_distances
    
    # d = euclidean_distances(firing_rate, firing_rate)
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

    e = int(np.floor(num))
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.scatter(np.log10(th[1:e]), np.log10(corr_integral[1:e]))

    print('distance calculation finish')






    """ Saving """
    from scipy.io import savemat
    matlab_save = {
        "t": t,
        "firing_rate":firing_rate,
        "bp_data": bp_data,
        "fr_isomap": fr_isomap,
        "h0": h0,
        "h1": h1,
        "h2": h2,
        "th": np.log10(th[1:e]),
        "corr_integral": np.log10(corr_integral[1:e])
    }

    savemat("..\\POSTPROCESSED_MAT_CODE\\DATA\\"+rat_numb+"_matlab_save.mat", matlab_save)

