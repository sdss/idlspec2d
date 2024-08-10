#!/usr/bin/env python3
import numpy as np

def find_nearest_indx(array, value):
    arr=False if isinstance(value, (int, float)) else True
       
    value = np.atleast_1d(value)
    indxs=np.zeros_like(value, dtype=int)
    array = np.asarray(array)
    for i, val in enumerate(value):
        indxs[i] = (np.abs(array - val)).argmin()
    if not arr:
        indxs = indxs[0]
    return indxs

