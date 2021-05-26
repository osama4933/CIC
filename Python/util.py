import math

def length(arr):
    return max(arr.shape)

def pol2cart(phi, rho):
    x = rho * math.cos(phi)
    y = rho * math.sin(phi)
    return(x, y)

import operator as op
from functools import reduce

def nchoosek(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom