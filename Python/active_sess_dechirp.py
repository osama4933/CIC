from param_configs import param_configs
import numpy as np
from util import length
from sym_to_data_ang import sym_to_data_ang
import math
import time

def active_sess_dechirp(x_1):
    #ACTIVE_SESS_DECHIRP Summary of this function goes here
    #   Detailed explanation goes here
    
    SF = param_configs(1)
    BW = param_configs(2)
    Fs = param_configs(3)
    N = int(2**SF)
    upsampling_factor = int(Fs/BW)
    
    DC = np.conj(sym_to_data_ang([1],N))
    DC_fft = np.fft.fft(DC)
    DC_upsamp = np.fft.ifft(np.concatenate([DC_fft[:N//2+1], np.zeros((upsampling_factor-1)*N), DC_fft[N//2 + 1:N]]))
    
    peak_gain = []
    uplink_wind = []
    n = []
    p = []
    last_wind = 0
    win_jump_factor = 3
    front_buf = 6*win_jump_factor
    back_buf = 3*win_jump_factor
    win_jump = math.floor(N*upsampling_factor/win_jump_factor)
    mov_thresh_wind = 1000*win_jump_factor
    mov_thresh = 0
    mov_thresh_rec  = []
    # t = time.perf_counter()
    # s1 = 0
    # s2 = 0
    # s3 = 0
    # s4 = 0
    for i in range(math.floor(length(x_1)/win_jump) - win_jump_factor):#length(DC))
        
        wind = x_1[i*win_jump : i*win_jump+(N*upsampling_factor)]
        # t1 = time.perf_counter()
        wind_fft = np.abs(np.fft.fft(wind * DC_upsamp)) # 17 sec
        wind_fft = wind_fft[np.concatenate([np.arange(N//2, dtype=int), \
                                            np.arange((N//2 + (upsampling_factor-1)*N), (upsampling_factor)*N, dtype=int)])]
        noise_floor = np.mean(wind_fft) # 4 sec
        # s1 += time.perf_counter() - t1
        n.append(noise_floor)
        fft_peak = wind_fft.max(0)
        
        # t1 = time.perf_counter()
        p.append(fft_peak)
        peak_gain.append(10*math.log10(fft_peak/noise_floor))
        
        if i+1 > mov_thresh_wind:
            # t2 = time.perf_counter()
            mov_thresh = 1.3*np.mean(peak_gain[-mov_thresh_wind + 1:])
            # s3 += time.perf_counter() - t2
            if mov_thresh > 6:
                mov_thresh = 6
        else:
            mov_thresh = 1.3*np.mean(peak_gain)
            if mov_thresh > 6:
                mov_thresh = 6
        # s2 += time.perf_counter() - t1
        mov_thresh_rec.append(mov_thresh)
        
        if peak_gain[-1] >= mov_thresh:
            if i+1 > last_wind:
                if i-back_buf < 1:
                    uplink_wind.append([1, i+1+front_buf])
                else:
                    uplink_wind.append([i+1-back_buf, i+1+front_buf])
                last_wind = uplink_wind[-1][1]
            elif i+1 <= last_wind:
                uplink_wind[-1][1] = i + 1 + front_buf
                last_wind = uplink_wind[-1][1]
    uplink_wind = np.array(uplink_wind)
    uplink_wind = uplink_wind[((uplink_wind[:,1]-uplink_wind[:,0]) != (front_buf + back_buf)).nonzero()[0],:]
    uplink_wind = uplink_wind[((uplink_wind[:,1]-uplink_wind[:,0]) != (front_buf + back_buf + 1)).nonzero()[0],:]
    uplink_wind = uplink_wind[((uplink_wind[:,1]-uplink_wind[:,0]) != (front_buf + back_buf - 1)).nonzero()[0],:]
    temp_link = uplink_wind
    uplink_wind = uplink_wind * win_jump
    
    if uplink_wind[-1,1] > length(x_1):
        uplink_wind[-1,1] = length(x_1)

    return uplink_wind