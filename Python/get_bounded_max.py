import numpy as np

def get_bounded_max(arr,up_thresh,low_thresh):
    (pnts_up,) = (arr < up_thresh).nonzero()
    (pnts_low,) = (arr > low_thresh).nonzero()
    pnts = np.intersect1d(pnts_up, pnts_low)
    return pnts
