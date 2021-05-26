import math
import numpy as np
from util import pol2cart

def sym_to_data_upsampled(symbols,N,Fs,BW):
    data = []
    accumulator = 0
    f0 = -BW/2
    fn = BW/2
    
    F = np.arange(f0 - 1, fn, BW/(N*(Fs/BW)))
    F = F[:-1]
    phase = 0
    for j in symbols:
        temp = []
        c = 1
        for i in range(N*(int(Fs)//int(BW))):
            phase = phase + ((2*math.pi*F[i])/Fs)
            polar_radius = 1

            [x, y] = pol2cart(phase, polar_radius)

            temp.append(complex(x, y))
            c = c+1
        data += temp
    return np.array(data)

