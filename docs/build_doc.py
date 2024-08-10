#!/usr/bin/env python3

import glob
import subprocess
import os.path as ptt
from os import environ, getenv
import argparse
import sys
import time

environ['DATABASE_PROFILE'] = 'READTHEDOCS'
if getenv('IDLUTILS_DIR') is None:
    environ['IDLUTILS_DIR'] = ptt.join(getenv('READTHEDOCS_VIRTUALENV_PATH'),'idlutils')


try:
    from pkg_resources import resource_filename
    bindir = resource_filename('boss_drp','../../bin/')
    prodir = resource_filename('boss_drp','../../pro/')
    docdir = resource_filename('boss_drp','../../docs/sphinx/')
except:
    file_path = ptt.realpath(__file__)
    print(file_path)
    bindir = ptt.join(ptt.dirname(ptt.dirname(file_path)),'bin/')
    prodir = ptt.join(ptt.dirname(ptt.dirname(file_path)),'pro/')
    docdir = ptt.join(ptt.dirname(ptt.dirname(file_path)),'docs/sphinx/')

mask = '\n.. _{name}:\n\n{name}\n{fmt}\n::\n \n    {doc}\n'


def headline(text, adorn='='):
    return text + '\n' + adorn*len(text)

def sec(out, sec_hdr, data, typestr):
    out.write(headline(sec_hdr, adorn='-')+'\n\n')
        
    out.write('.. contents::\n')
    out.write('    :depth: 3'+'\n')
    out.write('    :local:\n')
    out.write('    :class: this-will-duplicate-information-and-it-is-still-useful-here\n')
    out.write('    :backlinks: none\n\n')
        
    for cmd in data:
        cmd['doc'] = cmd['doc'].replace('\n','\n    ')
        out.write(mask.format(type='typestr', fmt = '^'*len(cmd['name']), **cmd))
    out.write('\n')
    return(out)

def build_docs():
    def filter(test,docstr):
        if test in docstr:
            dss = []
            for ds in docstr.split('\n'):
                if test not in ds: dss.append(ds)
            docstr = '\n'.join(dss)
        return(docstr)

    docs = {}

    docs['cmd'] = []
    for command in sorted(glob.glob(bindir+'/*'), key=ptt.basename):
        docstr = subprocess.getoutput(command+' -h')

        docstr = filter('Overriding default configuration',docstr)
        docstr = filter('PyFITSDeprecationWarning',docstr)
        docstr = filter('PyFITS is deprecated', docstr)
        docstr = filter('pyautogui does not seem to be available',docstr)
        docstr = filter('esutil not available!',docstr)
        docstr = filter('No slurm package installed:',docstr)
        docstr = filter('ERROR: dustmaps is not installed',docstr)
        docstr = filter('Environmental Varable IDLUTILS_DIR must be set',docstr)
        docstr = filter('WARNING: No SDSSDB access',docstr)
        docstr = filter('ERROR: No SDSSDB access',docstr)
        docstr = filter('No slurm package',docstr)
        docstr = filter('no gaiaxpy...!',docstr)
        docs['cmd'].append({'name':ptt.basename(command), 'doc': docstr})

    docs['idl'] = []
    for command in ['spreduce2d.pro','rm_combine_script.pro',
                    'spreduce1d_empca.pro','spcalib_qa.pro',
                    'spspec_target_merge.pro']:
        pf = glob.glob(prodir+'/*/'+command)
        if len(pf) == 0: continue
        docstr = []
        with open(pf[0],'r') as prof:
            docstr = prof.read()
        dss = []
        for ds in docstr.split('\n'):
            if ';------------------' in ds:
                break
            if len(ds.strip()) == 0:
                continue
            dss.append(ds)
        docstr = '\n'.join(dss)
        docs['idl'].append({'name':command, 'doc': docstr})

    with open(docdir+'doc.rst', 'w') as out:
        out.write(':tocdepth: 2\n\n')
        out.write('.. highlight:: none\n\n')
        out.write(headline('Full Command Documention') + '\n')
        out.write('Documented below are the primary commands used to run the BOSS Data Reduction Pipeline. However, there are numerous other routines included in this package, which are called by these commands and have their own internal documentation.\n\n')
        sec_hdr ='Full Bash and Python Command Usage'
        out = sec(out, sec_hdr, docs['cmd'], 'bin')

        sec_hdr ='IDL Command Usage'
        out = sec(out, sec_hdr, docs['idl'], 'idl')

        out.write('\n.. highlight:: defaults\n\n')

        out.write('\n.. End of document\n')

if __name__ == '__main__' :
    """
    Build BOSS DRP Documention
    """
    build_docs()
