
from boss_drp.field import field_to_string

import numpy as np

def fieldgroup(field, custom=False):
    field = str(field)
    if custom:
        if len(field.split('_')) == 1:
            return(field)
        else:
            return('_'.join(field.split('_')[:-1]))
    elif field.isnumeric():
        return(str(np.floor(int(field)/1000).astype(int)).zfill(3)+'XXX')
    elif field == '*':
        zfield = field_to_string(0)
        return(fieldgroup(zfield).replace('0','?'))
    else:
        return(field)
