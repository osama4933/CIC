import numpy as np
import math
from util import pol2cart

def sym_to_data_ang(symbols,N):
    data = []
    accumulator = 0
    
    for j in symbols:
        phase = -math.pi + ((j-1)*(2*math.pi/(N)))
        temp = []
        for i in range(N):
            accumulator = accumulator + phase
            polar_radius = 1
            [x, y] = pol2cart(accumulator, polar_radius)
            temp.append(x + 1j*y)
            phase = phase + (2*math.pi/(N))
        data += temp
    return np.array(data)


