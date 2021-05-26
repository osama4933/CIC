# NOTE: Up to date???
import numpy as np
from util import length
from get_max import get_max
from param_configs import param_configs
from sym_to_data_ang import sym_to_data_ang
import math

def UC_location_corr_DC_based(Data,DC_ind):
    # chirp variables
    SF = param_configs(1)
    BW = param_configs(2)
    Fs = param_configs(3)
    N = 2**SF
    upsampling_factor = Fs/BW

    # LORA pkt variables
    num_preamble = param_configs(4)
    num_sync = param_configs(5)
    num_DC = param_configs(6)
    num_data_sym = param_configs(7)
    # thresholds
    corr_threshold = param_configs(10)
    pnts_threshold = param_configs(11)

    DC = np.conj(sym_to_data_ang([1],N))

    ####################################
    if(len(DC_ind) == 0):
        return
    ## Find list of Potential Preambles from list of Downchirps Detected
    pot_pream_ind = []
    for i in range(len(DC_ind)):
        if DC_ind[i,0] - ((num_preamble+num_sync)*N) < 1:
            continue
        pot_pream_ind.append(np.arange(DC_ind[i,0] - ((num_preamble + num_sync)*N), DC_ind[i,0]- ((num_sync)*N) + 1, N))
    pot_pream_ind = np.array(pot_pream_ind)

    Upchirp_ind = []

    #  Cross Correlation with a Single UpChirp
    temp_wind = []
    for j in range(pot_pream_ind.shape[0]):
        if pot_pream_ind[j,0] - N <= 0:
            continue
        Data_buffer = []
        Data_buffer = Data[int(pot_pream_ind[j,0] - N) : int(pot_pream_ind[j,-1] + N)]
        temp = [0+0j]
        for i in range(length(Data_buffer) - length(DC)):
            temp.append(np.sum(np.multiply(Data_buffer[i + 1 : i + N + 1], DC[:N])) \
                / math.sqrt(np.sum(Data_buffer[i + 1: i + N + 1] * Data_buffer[i + 1 : i + N + 1].conj() ) * \
                np.sum( DC[:N] * DC[:N].conj())))
        temp_wind.append(temp)
    temp_wind = np.array(temp_wind)

    array_stack = []

    # iterate over each Downchirp Detected
    for m in range(temp_wind.shape[0]):
        
        n_samp_array = []
        peak_ind_prev = np.array([])
        for i in range(math.floor(length(temp_wind)/N)):
            # windowing Cross-Correlation Arrays correponsing to each pkt (window length N samples)
            wind = abs(temp_wind[m,i*N + 1 : (i+1) * N])
            peak_ind_curr = get_max(wind,corr_threshold,pnts_threshold)

            if length(peak_ind_prev) != 0 and length(peak_ind_curr) != 0:

                for j in range(length(peak_ind_curr)):
                    for k in range(length(peak_ind_prev)):
                        # check if combination of any two peaks in consecutive window are N samples apart
                        if abs(peak_ind_curr[j]) == abs(peak_ind_prev[k]):
                            n_samp_array.append(peak_ind_prev[k]+((i-1)*N)+(pot_pream_ind[m,0]-N-1) + 3) # add 3 to account for zero indexing
                        # This extracts a list of all peaks that are N samples apart
            peak_ind_prev = peak_ind_curr
        array_stack.append(n_samp_array)
    array_stack = np.array(array_stack)

    for m in range(len(array_stack)):
        n_samp_array = np.array(array_stack[m])
        
        for i in range(length(n_samp_array)):
            c = 0
            ind_arr = np.arange(n_samp_array[i] + N, n_samp_array[i] + N + (num_preamble-2)*N + 1, N)
            for j in range(len(ind_arr)):
                c = c + np.sum( n_samp_array == ind_arr[j] )
            # Find from the list all the peaks that appear consecutively for
            # more than 6 windows (Upchirp should give 8 peaks, N sampled apart)
            if c >= 6:
                if len(Upchirp_ind) != 0:
                    if np.sum(np.array(Upchirp_ind)[:,0] == n_samp_array[i]) != 1:
                        Upchirp_ind.append(np.concatenate([[n_samp_array[i]],ind_arr]))
                else:
                    Upchirp_ind.append(np.concatenate([[n_samp_array[i]], ind_arr]))
    # filter Upchirps that are with in 5 samples (same pkt detected multiple times due to peak energy spread)
    temp = []
    Upchirp_ind = np.array(Upchirp_ind)
    indices = np.concatenate([np.zeros((1,num_preamble)), Upchirp_ind])
    for i in range(1, indices.shape[0]):
        if len(temp) == 0:
            temp.append(indices[i,:])
        else:
            if min(abs(indices[i][0] - np.array(temp)[:,0])) > 5:
                temp.append(indices[i,:])
    temp = np.array(temp)
    Upchirp_ind = temp

    return [Upchirp_ind]