#!/usr/bin/env python3
import numpy as np

def find_nearest_indx(array, value):
    is_scalar = np.isscalar(value)  # check before converting
    array = np.asarray(array)
    value = np.atleast_1d(value)
    
    # Efficiently compute indices of closest values
    indxs = np.abs(array[:, np.newaxis] - value).argmin(axis=0)
    
    return indxs.item() if is_scalar else indxs

