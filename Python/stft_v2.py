import numpy as np
from util import length
import math
def stft_v2(Buffer, f_n):
    w_n = f_n
    Buff_len = length(Buffer) # Signal length
    # Frequency axis
    f_n = math.ceil(f_n/2) * 2+1
    Lf = (f_n - 1)/2
    # Time axis
    w_n = math.ceil(w_n/2) * 2+1
    Lw = (w_n - 1)/2
    # Initialize Spectrum to zero with appropriate size
    Spec = np.zeros((f_n,Buff_len), dtype=np.complex64) 
    ## Sliding window over signal
    for iterr in range(Buff_len):
        i_l = min([iterr, Lw, Lf])
        i_r = min([Buff_len-iterr - 1, Lw, Lf])
        iter_ind = np.arange(-i_l, i_r+1)
        ind1 = (iter_ind + iterr).astype('int')   # Time Indexing of the original signal
        ind = (iter_ind + Lf).astype('int')     # Frequency Indexing of the martix
        Spec[ind, iterr] = Buffer[ind1]
    # Computing FFT of stacked Windows
    Spec = np.fft.fft(Spec.T).T
    Spec = Spec*2/f_n # normalizing the FFTs
    return Spec