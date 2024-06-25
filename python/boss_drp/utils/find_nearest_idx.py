#!/usr/bin/env python3
import numpy as np

def find_nearest_indx(array, value):
    indxs=np.zeros_like(value)
    array = np.asarray(array)
    for i, val in enumerate(value):
        indxs[i] = (np.abs(array - val)).argmin()
    return indxs

