import numpy as np
import math
from get_bounded_max import get_bounded_max
from sym_to_data_ang import sym_to_data_ang
from param_configs import param_configs
from stft_v2 import stft_v2

def filter_false_postives(Data_stack,Upchirp_ind,Peak,FFO):
    # load parameters
    SF = param_configs(1)
    BW = param_configs(2)
    Fs = param_configs(3)
    N = int(2**SF)
    num_preamble = param_configs(4)
    num_sync = param_configs(5)
    num_DC = param_configs(6)
    S1 = param_configs(12)
    S2 = param_configs(13)

    DC = np.conj(sym_to_data_ang(np.ones(num_preamble),N))

    ffo = []
    Preamble_ind = []
    bin_offsets = []
    Data_out = []
    Peak_amp = []
    pream_peak_ind = []
    # row_ind contains the Frequency Bins around bin 1 where a
    # Preamble Peak can lie, assumption (-6*BW/2^SF <= Freq_off <= 6*BW/2^SF)
    row_ind = np.concatenate([range(N-5, N+1), range(7)])
    for m in range(Upchirp_ind.shape[0]):
        # extract 8 Preamble long window
        data_wind = Data_stack[m,int(Upchirp_ind[m,0]) - 1 : int(Upchirp_ind[m,0] + ((num_preamble )*N)) - 1]
        # Compute STFT to accurately find the Preamble Frequency Bin and any
        # bin offset (if any left)
        Spec = stft_v2(data_wind * DC,N)
        temp = []
        count = 0
        for i in np.nditer(row_ind):
            temp.append(np.sum(np.abs(Spec[i,:])))
            count = count + 1
        ind = np.argmax(temp)
        pream_peak_ind.append(row_ind[ind])
        sync1_ind = np.mod(pream_peak_ind[m] + 8,N)
        sync2_ind = np.mod(pream_peak_ind[m] + 16,N)
        if(sync1_ind == 0):
            sync1_ind = N
        if(sync2_ind == 0):
            sync2_ind = N
        sync_wind = Data_stack[m,int(Upchirp_ind[m,num_preamble-1] + N - 1): int(Upchirp_ind[m,num_preamble-1] + N + (num_sync*N) - 1)]

        sync_threshold_up = Peak[m,0] + 0.5*Peak[m,0]
        sync_threshold_low = Peak[m,0] - 0.5*Peak[m,0]
        
        sync_word1 = abs(np.fft.fft(sync_wind[:N] * DC[:N]))
        sync_word2 = abs(np.fft.fft(sync_wind[N:] * DC[:N]))
        if(sync_threshold_low < (2*np.sum(sync_word1)/N)):
            sync_threshold_low = (2*np.sum(sync_word1)/N)
        elif( sync_threshold_low < (2*np.sum(sync_word2)/N)):
            sync_threshold_low = (2*np.sum(sync_word2)/N)
        syn1_pnts = get_bounded_max(sync_word1,sync_threshold_up,sync_threshold_low)
        syn2_pnts = get_bounded_max(sync_word2,sync_threshold_up,sync_threshold_low)


        if(np.sum(syn1_pnts == sync1_ind) and np.sum(syn2_pnts == sync2_ind)):
            Preamble_ind +=  [Upchirp_ind[m,:]]
            if(pream_peak_ind[m] < N/2):
                bin_offsets.append(1 + (-np.mod(pream_peak_ind[m]+1,N))) # add one since p_p_i consists of indices
            else:
                print('second branch')
                bin_offsets.append(np.mod(N+2 - pream_peak_ind[m],N))
            Data_out += [Data_stack[m,:]]
            Peak_amp += [Peak[m,:]]
            ffo.append(FFO[m])
    Preamble_ind = np.array(Preamble_ind)
    Data_out = np.array(Data_out)
    Peak_amp = np.array(Peak_amp)
    return [Preamble_ind, bin_offsets, Data_out, Peak_amp, ffo]

