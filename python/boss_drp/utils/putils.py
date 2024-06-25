#!/usr/bin/env python3

import sys
import os
import subprocess
import shlex
import gzip
import time
"""
putils is a set of miscellaneous python tools.

Originally Written by Gary Kushner (LBL).  Nov 2009.  Latest update April 2010.
"""

def searchPath(name, paths):
    """Search a path for a name (file, direcory, link, etc).  Return the absolute
       path to the found file or None"""
    for path in paths:
        if os.path.exists(os.path.join(path, name)):
            return os.path.abspath(os.path.join(path, name))
    return None


def runCommand(cmd, echo=False, logCmd=None, prefix="", shell=False, limit=None, timeout=None):
    """Run a command with the option to asynchronously display or log output.
    
       If shell=False, the cmd needs to be a list, but if you pass in a string
       it will be parsed into a list.
    
       echo will echo output to stdout.
    
       logCmd is a function pointer to use to put the output into a log.
    
       Returns (return code, output)."""
    #prefix=prefix.encode()
    output = ""#.encode()
    
    #    Handle the command parsing
    if isinstance(cmd, str) and not shell:
        cmd = [c for c in shlex.split(cmd)]

    #    Call the process
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.STDOUT,
                         shell=shell)
    start = time.time()
    last = None
    i = 0
    #    Process output until process dies
    while True:
        l = p.stdout.readline().decode("utf-8")
        if not l: break
        output += l
        l = l[:-1]    # yea, only safe on unix...
        if limit is not None:
            if l == last:
                i +=1
                if i >= limit:
                    last = l
                    continue
            else:
                if i >= limit:
                    if echo:
                        print(prefix + f'..... {i-limit} duplicate log outputs removed .....')
                        print(prefix + last)
                    if logCmd != None:
                        logCmd(prefix + f'..... {i-limit} duplicate log outputs removed .....')
                        logCmd(prefix + last)
                i = 0
        if echo:
            print(prefix + l)
        if logCmd != None:
            logCmd(prefix + l)
        if limit is not None:
            last = l
        if timeout is not None:
            if time.time() - start > timeout:
                print(prefix + 'Timeout Command')
                if logCmd != None:
                    logCmd(prefix + 'Timeout Command')
                p.kill()
                break
    
    return (p.wait(), output)
    
def openRead(filename, mode = "r"):
    """Open a gzip or normal file for text reading.  Valid modes are 'r' and 'rb'"""
    
    gzSig = b'\x1f\x8b'
    
    if mode != 'r' and mode != 'rb':
        raise ValueError("Illegal mode: " + mode)
    
    f = open(filename, mode)
    
    try:
        if (f.read(2) == gzSig):
            f = gzip.open(filename, mode)
    finally:
        f.seek(0)
        
    return f
        

