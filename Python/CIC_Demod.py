import numpy as np
from get_bounded_max import get_bounded_max
from util import length, nchoosek
from get_max import get_max
from param_configs import param_configs
from sym_to_data_ang import sym_to_data_ang
import math

def CIC_Demod(Pream_ind,Rx_Buffer,Pream_ind_stack,Peak_amp,m):
    #DEMOD
    # CIC Demodulation
    # chirp variables
    SF = param_configs(1)
    N = int(2**SF)

    # LORA pkt variables
    num_preamble = param_configs(4)
    num_sync = param_configs(5)
    num_DC = param_configs(6)
    num_data_sym = param_configs(7)

    DC = np.conj(sym_to_data_ang([1], N))

    ###########################################################################
    Pream_ind_stack = np.append(Pream_ind_stack, Pream_ind_stack[:,num_preamble-1] + N)
    ###########################################################################
    
    # for each Preamble in the Pream_ind_stack, compute exact start and end
    # indices for all the data symbols
    frm_ind = []
    for i in range(Pream_ind_stack.shape[0]):
        frm_st = Pream_ind_stack[i] + (num_preamble*N) + (num_DC*N) + (num_sync*N)
        frm_ind.append([np.arange(frm_st, frm_st+((num_data_sym-1)*N)+1, N), \
                        np.arange(frm_st+N-1, frm_st+((num_data_sym)*N)+1, N)])
    frm_ind = np.array(frm_ind) - 1 # subtract 1 for index conversion

    # for the pkt to be demodulated, find indices for each data symbol
    data_ind = []
    en_arr = []
    sym_peak = []
    Data_frame_start = Pream_ind[0] + (num_preamble*N) + (num_DC*N) + (num_sync*N)
    Data_frame_end = Data_frame_start + (num_data_sym*N)
    frame_indices = np.array([(np.arange(Data_frame_start, Data_frame_start+((num_data_sym-1)*N)+1, N)),
                                    (np.arange(Data_frame_start+N-1, Data_frame_start+((num_data_sym)*N)+1, N))]).T.astype('int') - 1 # subtract 1 for zero indexing
    
    data_ind = []
    sym_peak = []
    symbols = []
    for k in range(num_data_sym):
        ## Find interfering Symbol Boundaries
        # for current demodulation window, find the chunks of interfering
        # symbols due to collisions by determining index overlaps
        ind = []
        sym_bnd = []
        for i in range(len(frm_ind)):
            if i != m:
                st = frm_ind[i][0]
                ed = frm_ind[i][1]
                newst = st[np.intersect1d((st > frame_indices[k,0]).nonzero() , (st < frame_indices[k,1]).nonzero())]
                newed = ed[np.intersect1d((ed > frame_indices[k,0]).nonzero() , (ed < frame_indices[k,1]).nonzero())]
                if len(newst) != 0:
                    sym_bnd.append(newst)
                if len(newed) != 0:
                    sym_bnd.append(newed)
        sym_bnd = np.array(sym_bnd)
        ## CIC Filtering
        if frame_indices[k, 1] >= len(Rx_Buffer):
            symbols.append(-3)
            continue
        # standard LoRa dechirping
        data_wind = Rx_Buffer[frame_indices[k,0]:frame_indices[k,1]+1] * DC
        data_fft = np.abs(np.fft.fft(data_wind))

        # scale the dmeodulation window with appropriate Gaussian to
        # suppress the interfering symbols at window ends
        sigma = 1
        WinFun = np.exp(-(1/(2*(sigma**2)))* np.linspace(-1,1,N) ** 2)
        WinFun = WinFun / (math.sqrt(2*math.pi)*sigma)
        temp_wind = data_wind * WinFun
        sym_bnd = np.mod(sym_bnd - frame_indices[k,0],N)
        intf_wind = []
        # n_pnt is the fft Factor - fft(Signal, n_pnt * length(Signal))
        n_pnt = 4
        for i in range(len(sym_bnd)):
            buff = np.zeros((2,n_pnt*N), dtype=np.complex64)
            sym_bnd_i = int(sym_bnd[np.unravel_index(i, sym_bnd.shape, 'F')])
            buff[0,:sym_bnd_i - 1] = temp_wind[:sym_bnd_i - 1]
            buff[0,:] = np.abs(np.fft.fft(buff[0,:], n_pnt*N)) / np.sqrt(np.sum(np.abs(buff[0,:]) ** 2)) # div OKAY
            buff[1,sym_bnd_i-1:N] = temp_wind[sym_bnd_i-1:N]
            buff[1,:] = np.abs(np.fft.fft(buff[1,:],n_pnt*N)) / np.sqrt(np.sum(np.abs(buff[1,:]) ** 2))
            intf_wind = buff
        # CIC's min operation to suppress interfering symbol's Peaks and
        # then finding candidate symbols
        intf_wind = np.array(intf_wind)
        intf_wind_min_fft = np.amin(intf_wind,axis=0) if len(intf_wind) != 0 else []
        pot_sym_cic = get_max(intf_wind_min_fft,4*np.sum(intf_wind_min_fft)/(n_pnt*N),n_pnt*N)
        pot_sym_cic = np.ceil(pot_sym_cic/n_pnt)
        # Out of all peaks, find ones with in range of +- 0.5*Preamble Peak
        PwrFctr = 0.5
        PwrFlr = 4 # may try 3
        up_thresh = (Peak_amp[0] + PwrFctr*Peak_amp[0])
        low_thresh = (Peak_amp[0] - PwrFctr*Peak_amp[0])
        if(low_thresh < (PwrFlr*np.sum(data_fft)/N)): #1
            low_thresh = (PwrFlr*np.sum(data_fft)/N)
        pot_sym_pf = get_bounded_max(data_fft,up_thresh,low_thresh)
        ## Filtering Preamble of interfering Packets
        # Filter out Peaks with in current window that are appearinf repeatedly
        # in 3 consecutive windows
        if not (frame_indices[k,1] + N > len(Rx_Buffer) or frame_indices[k,1] + 2*N > len(Rx_Buffer)):
            data_wind_next_1 = Rx_Buffer[frame_indices[k,0] + N : frame_indices[k,1] + N+1] * DC
            data_wind_prev_1 = Rx_Buffer[frame_indices[k,0] - N : frame_indices[k,1] - N+1] * DC
            data_wind_next_2 = Rx_Buffer[frame_indices[k,0] + 2*N : frame_indices[k,1] + 2*N+1] * DC
            data_wind_prev_2 = Rx_Buffer[frame_indices[k,0] - 2*N : frame_indices[k,1] - 2*N+1] * DC
            temp_next_1 = np.abs(np.fft.fft(data_wind_next_1,N))
            temp_prev_1 = np.abs(np.fft.fft(data_wind_prev_1,N))
            temp_next_2 = np.abs(np.fft.fft(data_wind_next_2,N))
            temp_prev_2 = np.abs(np.fft.fft(data_wind_prev_2,N))    
            next_wind_sym_1 = get_bounded_max(temp_next_1,4*np.sum(temp_next_1)/N,N) # TODO changed these four lines
            next_wind_sym_2 = get_bounded_max(temp_next_2,4*np.sum(temp_next_2)/N,N)
            prev_wind_sym_1 = get_bounded_max(temp_prev_1,4*np.sum(temp_prev_1)/N,N)
            prev_wind_sym_2 = get_bounded_max(temp_prev_2,4*np.sum(temp_prev_2)/N,N)

            temp = []
            for i in range(length(pot_sym_pf)):
                if( (np.sum(pot_sym_pf[i] == prev_wind_sym_1) and np.sum(pot_sym_pf[i] == next_wind_sym_1))\
                        or (np.sum(pot_sym_pf[i] == prev_wind_sym_2) and np.sum(pot_sym_pf[i] == prev_wind_sym_1))\
                        or (np.sum(pot_sym_pf[i] == next_wind_sym_1) and np.sum(pot_sym_pf[i] == next_wind_sym_2)) ):
                    pass                
                else:
                    temp.append(pot_sym_pf[i])
            pot_sym_pf = np.array(temp) + 1 # convert from indices to data
        ##  Freq. Offset Filtering
        # since we have removed Frequency Offset from pkt under consideration 'Pream_ind'
        # and chosen the right downsampled buffer, the true symbol peak should
        # be the most crisp one (either use this or Choir Module(next), results should be almost similar)
        # and the interfering symbols peak may or may not be crisp
        temp = []
        for i in range(length(pot_sym_pf)):
            if sum(pot_sym_pf[i] + 1 == pot_sym_pf) or sum(pot_sym_pf[i] - 1 == pot_sym_pf): #TODO what is going on here???
                pass
            else:
                temp.append(pot_sym_pf[i])
        pot_sym = np.array(temp)
        ##  Choir Module
        # npnt = 16
        # data_fft_npnt = abs(np.fft.fft(data_wind,npnt*N))
        # FO_thresh = 0.25
        # sym_FO = []
        # temp = []
        # for i in range(length(pot_sym_pf)):
        #     ind = []
        #     if(pot_sym_pf[i] == 1):
        #         ind = np.concatenate([np.arange((N*npnt) - (npnt/2) + 1, (N*npnt)), np.arange((((pot_sym_pf[i]-1) * npnt) + 1) + (npnt/2), N*npnt)])
        #     else:
        #         ind = np.arange((((pot_sym_pf[i]-1) * npnt) + 1) - (npnt/2), (((pot_sym_pf[i]-1) * npnt) + 1) + (npnt/2) + 1)
        #     # ind OKAY
        #     ind = ind.astype('int') - 1 # convert back from data to indices
        #     a = np.argmax(data_fft_npnt[ind])
        #     sym_FO.append(abs(a - ((npnt/2)+1))/npnt)
        #     if(sym_FO[-1] < FO_thresh):
        #         temp.append(pot_sym_pf[i])
        # pot_sym = np.array(temp)
        #################################################################
        ## Make the final decision
        b = []
        if(length(sym_bnd) == 0):
            # if there is no symbol colliding with current demod window
            if(len(pot_sym) == 0):
                symbols.append(np.argmax(data_fft) + 1)
            else:
                # choose peak closest in height to Preamble Peak
                dist = abs(data_fft[pot_sym - 1] - (up_thresh + low_thresh)/2)
                b = np.argmin(dist)
                symbols.append(pot_sym[b])
        else:
            # if symbols are colliding with current demod window
            fin_sym = (np.intersect1d(pot_sym_cic,pot_sym)).astype('int')
            ##  Final Decision
            if(length(fin_sym) == 0):
                if length(pot_sym_cic) == 0 and length(pot_sym) != 0:
                    # choose peak closest in height to Preamble Peak
                    dist = abs(data_fft[pot_sym - 1] - (up_thresh + low_thresh)/2)
                    b = np.argmin(dist)
                    symbols.append(pot_sym[b])
                    
                elif length(pot_sym) == 0 and length(pot_sym_cic) != 0:
                    # make decision based on CIC's windows, correct symbol should
                    # have lowest std as it appears in all windows
                    sdev = np.std(intf_wind[:,(n_pnt * pot_sym_cic).astype('int')])
                    b = np.argmin(sdev)
                    symbols.append(pot_sym_cic[b])
                    
                elif length(pot_sym) == 0 and length(pot_sym_cic) == 0:
                    symbols.append(np.argmax(data_fft) + 1) # add one to account for zero-indexing
                else:
                    # choose peak closest in height to Preamble Peak
                    dist = abs(data_fft[pot_sym - 1] - (up_thresh + low_thresh)/2)
                    b = np.argmin(dist)
                    symbols.append(pot_sym[b])
                    
            else:
                # if intersection yields some candidates then decide based on partial STFT as following
                
                ##  Stft
                # find stft 2D matrix of follwoing dimensions
                # N frequency (rows)  x  [1 : avg_pnts      N - avg_pnts : N] (columns)
                # i.e. finding Spectrum of first 10 and last 10 time samples (Dont need to compute whole Spectrum)
                avg_pnts = 10 # number of start and end time samples to average over
                G_wind1 = data_wind
                Spec = np.zeros((N, 2*(avg_pnts+1)), dtype=np.complex64)
                for i in range(avg_pnts+1):
                    Spec[:int(N//2) + i,i] = G_wind1[:int(N//2) + i]
                    Spec[int(N//2) - (avg_pnts - i)-1:N,i+1+avg_pnts] = G_wind1[int(N//2) - (avg_pnts - i)-1:N]
                Spec = np.fft.fft(Spec, axis=0)
                # the amplitude difference at the start and end of the spectrum
                # for correct symbol should be minimum, (non-interfering symbol's
                # frequency track appears continuosly for all the columns)
                spec_slice1 = np.abs(Spec[fin_sym-1,:(avg_pnts+1)]) # TODO why minus 2??
                freq_amp = np.amin(spec_slice1,axis=1) if len(spec_slice1) != 0 else []
                spec_slice2 = np.abs(Spec[fin_sym-1,-avg_pnts-1:])
                freq_amp_end = np.amin(spec_slice2,axis=1) if len(spec_slice2) != 0 else []
                dif = np.abs(freq_amp - freq_amp_end)
                b = np.argmin(dif)
                symbols.append(fin_sym[b])
                
    return symbols