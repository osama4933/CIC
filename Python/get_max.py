import numpy as np

def get_max(arr,threshold,num_pnts):
    out = []
    if len(arr) == 0:
        return np.array(out)
    for i in range(num_pnts):
        [a, b] = arr.max(0), arr.argmax(0)
        if(a<threshold):
            return np.array(out)
        else:
            out.append(b)
            arr[b] = 0
    return np.array(out)

