from sym_to_data_ang import sym_to_data_ang
from sym_to_data_upsampled import sym_to_data_upsampled
from DC_location_correlation import DC_location_correlation
from UC_location_corr_DC_based import UC_location_corr_DC_based
from dnsamp_buff import dnsamp_buff
from filter_false_postives import filter_false_postives
from active_sess_dechirp import active_sess_dechirp
from param_configs import param_configs
from util import length
from CIC_Demod import CIC_Demod
import numpy as np
import math
import time

num_nodes = 15
num_data_sym = 28
in_file_path = '15tx_0.5lm_1' # all relative to main.py
symbols_out_path = 'symbols_out.txt'
peaks_out_path = 'peaks.txt'

## Loading variables
# chirp variables
SF = param_configs(1)
BW = param_configs(2)
Fs = param_configs(3)
N = int(2**SF)
upsampling_factor = int(Fs/BW)
Ts = 1/Fs

# LORA pkt variables
num_preamble = param_configs(4)
num_sync = param_configs(5)
num_DC = param_configs(6)
num_data_sym = param_configs(7)
preamble_sym = 1
pkt_len = num_preamble + num_sync + num_DC + num_data_sym
num_samples = pkt_len * N
# load true symbols
sym = np.loadtxt('symbols.txt') # TODO: Figure this out! Probably want this to be a .mat file

# Generating a Downchirp
DC = np.conj(sym_to_data_ang([1],N))

# parse complex data
x_1 = np.fromfile(in_file_path, dtype=np.complex64)

# file duration
t = np.arange(len(x_1)) / Fs

x_1 = x_1[:int(math.floor(length(x_1) / upsampling_factor) * upsampling_factor)]

x_1_dnsamp = x_1[::int(upsampling_factor)]        
file_dur = length(x_1) / Fs

uplink_wind = active_sess_dechirp(x_1) - 1 # subtract 1 for indexing
print('Detected ', len(uplink_wind), ' active sessions')

demod_sym_stack = []
Peaks = []
for m in range(uplink_wind.shape[0]):
    print('Packet ' + str(m+1) + '/' + str(uplink_wind.shape[0]))
    # DC correlations to find LoRa pkts out of collision
    temp_buff = []
    temp_buff = x_1[int(uplink_wind[m, 0]) : int(uplink_wind[m, 1]) + 1]
    temp_buff = temp_buff[:(len(temp_buff)//upsampling_factor)*upsampling_factor]
    DC_ind = DC_location_correlation(temp_buff[::upsampling_factor])
    print('Found ', len(DC_ind), ' Downchirps in current Active session')
    if (DC_ind.shape[0] == 0):
        continue

    ##      UC correlation to filter false positives and frequency offset & Packets' SNR estimation
    # All possible downsampled Buffers with different starting sample for downsampling
    Data_freq_off = np.array([])
    Rx_Buff_dnsamp = []
    for i in range(upsampling_factor):
        Rx_Buff_dnsamp.append(temp_buff[i::upsampling_factor])
    Rx_Buff_dnsamp = np.array(Rx_Buff_dnsamp)
    [Upchirp_ind] = UC_location_corr_DC_based(temp_buff[::upsampling_factor], DC_ind)
    if (Upchirp_ind.shape[0] == 0):
        continue

    # for each Preamble detected, Choosing the correct downsampled buffer with any frequency offset
    # been compensated and determining Preamble Peak heights to be used later for Power filtering
    [Data_freq_off, Peak, Upchirp_ind, FFO] = dnsamp_buff(Rx_Buff_dnsamp, Upchirp_ind)
    if (Upchirp_ind.shape[0] == 0):
        continue

    # Filter False Positives based on 2-SYNC Words detected 
    [Preamble_ind, bin_offsets, Data_out, Peak_amp, FFO] = filter_false_postives(Data_freq_off, Upchirp_ind, Peak, FFO)
    # filter preambles that are with in 5 samples (same pkt detected twice due to Correlation peak energy spread)
    temp = []
    temp_data = []
    temp_peaks = []
    indices = np.vstack([[np.zeros(num_preamble)], Preamble_ind])
    Data = np.vstack([[np.zeros(Data_out.shape[1])], Data_out])
    peaks = np.vstack([[np.zeros(Peak_amp.shape[1])], Peak_amp])
    Peak_amp = []
    for i in range(1, indices.shape[0]):
        if (len(temp) == 0):
            temp.append(indices[i, :])
            temp_data.append(Data[i, :])
            temp_peaks.append(peaks[i, :])
        else:
            if (min(abs(indices[np.unravel_index(i, indices.shape, 'F')] - np.array(temp)[:, 1])) > 10):
                temp.append(indices[i, :])
                temp_data.append(Data[i, :])
                temp_peaks.append(peaks[i, :])
    Pream_ind = np.array(temp)
    Data_out = np.array(temp_data)
    Peak_amp = np.array(temp_peaks)
    ##  Data Demodulation using CIC
    print('Found ', len(Pream_ind), ' Preambles in current Active session')
    demod_sym = []
    for j in range(Pream_ind.shape[0]):
        demod_sym.append(CIC_Demod(Pream_ind[j, :], Data_out[j, :], Pream_ind, Peak_amp[j, :], j))
        demod_sym[j] = np.mod(demod_sym[j] + bin_offsets[j] - 2, N)
    demod_sym_stack.append(demod_sym[0])
    Peaks.append(Peak_amp)

# Saving symbols and peaks to disk
demod_sym_stack = np.array(demod_sym_stack)
Peaks = np.array(Peaks[0])
np.savetxt(symbols_out_path, demod_sym_stack, fmt="%i")
np.savetxt(peaks_out_path, Peaks)
# Calculating SER for the file
ser = np.sum(np.tile(sym, (demod_sym_stack.shape[0], 1)) != demod_sym_stack) \
    / np.prod(demod_sym_stack.shape)
print('******************************************')
print('Symbol Error Rate for this File = ', ser)
print('*****************FINISHED*****************')