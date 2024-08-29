#!/usr/bin/env python3
from boss_drp.utils.lock import lock, unlock

import glob
import hashlib
import os
import os.path as ptt
from pathlib import Path
import time


def compute_sha1(file_path):
    """Computes the SHA-1 hash of the given file."""
    file = Path(file_path)
    sha1 = hashlib.sha1(file.read_bytes())
    return sha1.hexdigest()


def create_hash_line(file, path):
    filer= ptt.relpath(file, start=path)
    out = '{}  {}'.format(compute_sha1(file), filer)
#    out = '{}  {}'.format(hsh.hexdigest(), file.name)
    return out


def create_hash(path):
    path = ptt.abspath(path)
    output_file = ptt.join(path,'{}.sha1sum'.format(ptt.basename(path)))
    i = 0
    
    if lock(output_file, pause=5, niter=6):
        try:
            hash = []
            files = list(filter(ptt.isfile, glob.glob(ptt.join(path,'**'),recursive = True)))
            files = [x for x in files if '.sha1sum' not in x]
            files.sort(key=ptt.getmtime)
            for f in files:
                if ptt.isfile(f):
                    hash.append(create_hash_line(f, path)+'\n')
                    
            with open(output_file, 'w') as out:
                for h in hash:
                    out.write(h)
        finally:
            return not unlock(output_file)
    else:
        return(True)
    
    return(False)
    

def check_hash(data_dir, verbose=True):
    path = ptt.abspath(data_dir)
    output_file = ptt.join(path,'{}.sha1sum'.format(ptt.basename(path)))
    
    if not ptt.exists(output_file):
        if verbose:
            print('Missing Hash File: '+output_file)
        return(False)
    with open(output_file, 'r') as test:
        sh1 = test.readlines()
    sha1sums = {}
    valid = True
    for s1 in sh1:
        parts = s1.strip().split()
        if len(parts) == 2:
            sha1, filename = parts
            sha1sums[filename] = sha1
    for filename, expect_sha1 in sha1sums.items():
        filename = filename.replace('/',os.sep)
        if compute_sha1(ptt.join(data_dir, filename)) != expect_sha1:
            valid = False
            print(f"{filename}: FAILED (expected {expected_sha1}, got {computed_sha1})")
    if valid is True:
        print('All hashes are OK')
    return(valid)
    
    
