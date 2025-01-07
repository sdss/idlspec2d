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


def runCommand(cmd, echo=False, logCmd=None, prefix="", shell=False,
              limit=None, timeout=None, errCmd=None):
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

#    if errCmd is None:
#        errCmd = logCmd

    #    Call the process
    if errCmd is None:
        print('ets')
        err = subprocess.STDOUT
    else:
        subprocess.PIPE
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = err,
                         shell=shell)
    start = time.time()
    last = None
    i = 0
    #    Process output until process dies
    while True:
        # Read both stdout and stderr
        stdout_line = p.stdout.readline().decode("utf-8")
        try:
            stderr_line = p.stderr.readline().decode("utf-8")
        except:
            stderr_line = None

        # Break when both stdout and stderr are done
        if errCmd is not None:
            if not stdout_line and not stderr_line:
                break
        else:
            if not stdout_line:
                break

        # Capture stdout and stderr
        if stdout_line:
            output += stdout_line
        if stderr_line:
            output += stderr_line

        # Process the output
        if stdout_line:
            stdout_line = stdout_line.strip()  # Safe for Unix

            if limit is not None:
                if stdout_line == last:
                    i += 1
                    if i >= limit:
                        last = stdout_line
                        continue
                else:
                    if i >= limit:
                        if echo:
                            print(prefix + f'..... {i-limit} duplicate log outputs removed .....')
                            print(prefix + last)
                        if logCmd is not None:
                            logCmd(prefix + f'..... {i-limit} duplicate log outputs removed .....')
                            logCmd(prefix + last)
                    i = 0

            if echo:
                print(prefix + stdout_line)
            if logCmd is not None:
                logCmd(prefix + stdout_line)

            if limit is not None:
                last = stdout_line

        if stderr_line:
            stderr_line = stderr_line.strip()  # Safe for Unix

            if echo:
                print(prefix + stderr_line)
            if logCmd is not None:
                logCmd(prefix + stderr_line)

        # Timeout check
        if timeout is not None:
            if time.time() - start > timeout:
                print(prefix + 'Timeout Command')
                if logCmd is not None:
                    logCmd(prefix + 'Timeout Command')
                p.kill()
                break
    return p.wait(), output

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
        

