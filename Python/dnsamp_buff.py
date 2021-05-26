import numpy as np
from util import length
from stft_v1 import stft_v1
from param_configs import param_configs
from sym_to_data_ang import sym_to_data_ang
import math, cmath

def dnsamp_buff(Data_stack,Upchirp_ind):
    # load Parameters
    SF = param_configs(1)
    BW = param_configs(2)
    Fs = param_configs(3)
    N = int(2**SF)
    num_preamble = param_configs(4)
    num_sync = param_configs(5)
    num_DC = param_configs(6)
    DC = np.conj(sym_to_data_ang([1],N))

    ####################################
    ##  Compute and Correct Frequency Offsets for each Preamble Detected in each Data_stack and Find the Peak Statistics needed for demodulation

    Up_ind = []
    peak_amp = []
    Data_buff = []
    ffo = []
    FFO = []
    # n_pnt is the fft Factor - fft(Signal, n_pnt * length(Signal))
    n_pnt = 16
    peak_stats = []
    # iterate over all Upchirps that qualified 8 consecutive Peak condition
    for k in range(Upchirp_ind.shape[0]):
        if(Upchirp_ind[k,0] - N <= 0):
            peak_stats.append([])
            continue
        inn = []
        k_peak_stats = []
        Data_freq_off = []
        # iterate overall downsampled buffers
        for m in range(Data_stack.shape[0]):
            data_wind = []
            data_fft = []
            freq_off = []
            # ind_temp contains the Frequency Bins around bin 1 where a
            # Preamble Peak can lie
            ind_temp = np.concatenate([np.arange(5*n_pnt), np.arange((N*n_pnt)-(4*n_pnt)-1, (N*n_pnt))])
            # iterate over all Preambles
            c = []
            for j in range(num_preamble):
                data_wind = Data_stack[m,int(Upchirp_ind[k,0]) - 1 : int(Upchirp_ind[k,0] + (num_preamble*N) -1)]
                data_fft.append(abs(np.fft.fft(data_wind[((j)*N):((j+1)*N)] * DC[:N],n_pnt*N)))
                
                c.append(data_fft[j][ind_temp].argmax(0))
                c[j] = ind_temp[c[j]] + 1
                # Handle -ve and +ve Frequency Offsets Accordingly
                if(c[j] > (n_pnt*N)/2):
                    freq_off.append(( (N*n_pnt) - c[j] ) / n_pnt)
                else:
                    freq_off.append(-1*( c[j] - 1 ) / n_pnt)
            # average the frequency offset of 6 middle Preambles
            freq_off = np.sum( freq_off[1:7] ) / (num_preamble - 2)
            ffo.append(freq_off)
            # Correct for the Frequency Offset in corresponding Data_Stack
            Data_freq_off.append(Data_stack[m,:] * np.exp( (1j*2*math.pi*(freq_off / N)) * np.arange(1, length(Data_stack[m,:]) + 1) ))
            # ind_temp contains the Frequency Bins around bin 1 where a
            # Preamble Peak can lie, assumption (-5*BW/2^SF <= Freq_off <= 5*BW/2^SF)
            ind_temp = np.concatenate([range(5), range(N-4, N)])
            a = []
            c = []
            data_wind = []
            data_fft = []
            # for the frequency offset corrected Data Stack, find FFT of Preamble to get Peak Statistics 
            for j in range(num_preamble):
                data_wind = Data_freq_off[m][int(Upchirp_ind[k,0]) - 1 : int(Upchirp_ind[k,0] + (num_preamble*N)) - 1]
                data_fft.append(abs(np.fft.fft(data_wind[(j)*N : (j+1)*N] * DC[:N],N)))
                [aj,cj] = data_fft[j][ind_temp].max(0), data_fft[j][ind_temp].argmax(0)
                a.append(aj); c.append(cj)

                c[j] = ind_temp[c[j]]
            k_peak_stats.append([np.mean(a), np.var(a, ddof=1), np.std(a, ddof=1)])
            
            ##  Find the Right Data_stack to work with
            # first find the stft of given stack at the Preamble Region,
            # Spec is a 2D Matrix, rows - Freq. Bins & Col. - Time Samples
            Spec = stft_v1(Data_freq_off[m][int(Upchirp_ind[k,0] - N)-1:int(Upchirp_ind[k,-1] + N - 1 - N)],N,DC[:N],0,0)
            temp = []
            freq_track_qual = []
            pream_peak_ind = []
            adj_ind = []
            # row_ind contains the Frequency Rows around bin 1 where a
            # Preamble Peak can lie
            row_ind = np.concatenate([range(N-6,N), range(0,6)])
            count = 1
            for i in np.nditer(row_ind):
                temp.append(np.sum(np.abs(Spec[i,:])))
                count = count + 1
            temp = np.array(temp)
            # Frequency Track in row containing Preamble should have
            # maximum energy
            ind = temp.argmax(0)
            pream_peak_ind = row_ind[ind]
            # Find row indices for Preamble row + 1 & - 1
            adj_ind = np.array([np.mod(pream_peak_ind-1+1,N), np.mod(pream_peak_ind+1+1,N)]) # plus 1 for index conversion
            if(np.sum(adj_ind == 0) == 1):
                adj_ind[(adj_ind == 0).nonzero()] = N
            # A good quality frequency track for a preamble is one that has
            # least energy leakage in adjacent rows (this promises very sharp FFT peaks)
            adj_ind -= 1 # subtract 1 to convert back to Python indices
            freq_track_qual = ( np.sum(np.abs(Spec[pream_peak_ind,:])) - np.sum(np.abs(Spec[adj_ind[0],:])) ) + ( np.sum(np.abs(Spec[pream_peak_ind,:])) - np.sum(np.abs(Spec[adj_ind[1],:])) )
            inn.append(freq_track_qual)
        inn = np.array(inn)
        peak_stats.append(k_peak_stats)
        Data_freq_off = np.array(Data_freq_off)
        # choosing the best Data_stack based on maximum energy difference from
        # adjacent bins
        b = inn.argmax(0)
        # output frequency offset corrected buffer with relevant, Peak
        # statistics and frequency offsets
        Data_buff.append(Data_freq_off[b,:])
        FFO.append(ffo[b])
        peak_amp.append(peak_stats[k][b])
        Up_ind.append(Upchirp_ind[k,:])
    Data_buff = np.array(Data_buff)
    Up_ind = np.array(Up_ind)
    peak_amp = np.array(peak_amp)
    return [Data_buff, peak_amp, Up_ind, FFO]
