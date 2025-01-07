#!/usr/bin/env python3
import mmap


def grep(filepath, grepstr, line=False):
    with open(filepath, 'rb', 0) as f:
        s = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        if s.find(grepstr.encode('UTF-8')) != -1:
            if line:
                s.seek(0)  # Reset position to start of the file
                for ln in s.read().decode('UTF-8').splitlines():
                    if grepstr in ln:
                        return (True, ln)
            return True  # Found, but line=False
    return (False, None) if line else False
