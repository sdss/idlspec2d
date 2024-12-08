#!/usr/bin/env python3
import mmap

def grep(filepath, grepstr):
    with open(filepath, 'rb', 0) as f:
        s = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        if s.find(grepstr.encode('UTF-8')) != -1:
            return(True)
    return(False)
