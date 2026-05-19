#!/usr/bin/env python3
import re
import numpy as np

def match(array, value):
    if '*' in value:
        if r'[\w]*' not in value:
            value = value.replace('*',r'[\w]*')
    else:
        value = '('+value+'$)|('+value+r'\W)'
    r = re.compile(value, re.IGNORECASE)
    ret = np.full(len(array), False)
    idx = [i for i, x in enumerate(array) if r.search(x)]
    ret[idx] = True
    return(ret)

