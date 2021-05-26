import numpy as np
import math
from util import length
from get_max import get_max
from param_configs import param_configs
from sym_to_data_ang import sym_to_data_ang
import time

def DC_location_correlation(Rx_Buffer):
#   Detecting downchirp

    SF = param_configs(1)
    BW = param_configs(2)
    Fs = param_configs(3)
    N = int(2**SF)
    num_DC = param_configs(6)
    # thresholds
    corr_threshold = param_configs(8)      # Threshold above which we extract all Correlation peaks
    pnts_threshold = param_configs(9)      # Max. # of peaks to extract from Corrrelation Plot

    DC = np.conj(sym_to_data_ang([1],N))


    Downchirp_ind = []
    Cross_Corr = []

    st = time.perf_counter()
    ###################################################
    # Cross Correlation with a Single downchirp
    # OLD:
    for i in range(length(Rx_Buffer) - length(DC) - 1):
        Cross_Corr.append(np.sum(Rx_Buffer[ i : i + N ] * DC.conj()) \
                / math.sqrt(np.sum( Rx_Buffer[ i : i + N ] * Rx_Buffer[ i : i + (N) ].conj() ) *
                np.sum( DC * DC.conj())))
    Cross_Corr = np.array(Cross_Corr)
    Cross_Corr = Cross_Corr[np.isfinite(Cross_Corr)]
    corr_threshold =  4*np.sum(np.abs(Cross_Corr))/length(Cross_Corr)
    ###################################################
    n_samp_array = []
    peak_ind_prev = np.array([])
    for i in range(math.floor(length(Cross_Corr)/N)):
        # windowing Cross-Correlation (window length N samples)
        wind = np.abs(Cross_Corr[i*N : (i+1) * N])
        # Extract Multiple Correlation Peaks
        peak_ind_curr = get_max(wind,corr_threshold,pnts_threshold)
        if(length(peak_ind_prev) != 0 and length(peak_ind_curr) != 0):
            for j in range(length(peak_ind_curr)):
                for k in range(length(peak_ind_prev)):
                    # check if combination of any two peaks in consecutive window are N samples apart
                    if(peak_ind_curr[j] == peak_ind_prev[k]):
                        n_samp_array += [peak_ind_prev[k]+((i-1)*N) + 1, peak_ind_curr[j]+((i)*N) + 1]
                    # This extracts a list of all peaks that are N samples apart
        peak_ind_prev = peak_ind_curr

    n_samp_array = np.array(n_samp_array)
    for i in range(length(n_samp_array)):
        c = 0
        ind_arr = [n_samp_array[i], n_samp_array[i] + N]
        
        for j in range(len(ind_arr)):
            c = c + np.sum( n_samp_array == ind_arr[j] )
        # Find from the list all the peaks that appear consecutively for
        # more than 2 windows (Downchirp should give 3 peaks, N sampled apart)

        if( c >= 2 ):
            Downchirp_ind += [ind_arr]
    # filter Downchirps that are with in 3 samples (same pkt detected twice due to peak energy spread)
    # filter Downchirps that are with in 3 samples (same pkt detected twice due to peak energy spread)
    temp = []
    indices = [np.zeros(math.floor(num_DC))] + Downchirp_ind
    indices = np.array(indices)
    for i in range(1,indices.shape[0]):
        if(len(temp) == 0):
            temp.append(indices[i])
        else:
            if( np.min(np.abs(indices[np.unravel_index(i, indices.shape, 'F')] - np.array(temp)[:, 0])) > 3 ):
                temp.append(indices[i])
    Downchirp_ind = np.array(temp)
    return Downchirp_ind

