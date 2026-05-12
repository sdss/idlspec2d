#!/usr/bin/env python3

import glob
import subprocess
import os.path as ptt
from os import environ, getenv
import shutil
from pathlib import Path

environ['DATABASE_PROFILE'] = 'READTHEDOCS'
if getenv("IDLUTILS_DIR") is None:
    environ["IDLUTILS_DIR"] = str(
        Path(getenv("READTHEDOCS_VIRTUALENV_PATH")) / "idlutils"
    )

try:
    from pkg_resources import resource_filename
    bindir = resource_filename('boss_drp','../../bin/')
    prodir = resource_filename('boss_drp','../../pro/')
    docdir = resource_filename('boss_drp','../../docs/sphinx/')
except:
    file_path = Path(__file__).resolve()
    print(file_path)
    bindir = Path(file_path).parent.parent / "bin"
    prodir = Path(file_path).parent.parent / "pro"
    docdir = Path(file_path).parent.parent / "docs" / "sphinx"

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

def build_legacy_docs(py=True, bash=False, idl = False):
    if (not py) and (not bash) and (not idl):
        return 
    
    def filter(test,docstr):
        if test in docstr:
            dss = []
            for ds in docstr.split('\n'):
                if test not in ds: dss.append(ds)
            docstr = '\n'.join(dss)
        return(docstr)

    docs = {}
    docs['cmd'] = []
    for command in sorted(bindir.glob("*"), key=lambda p: p.name):
        if command.suffix.lower() == ".bash":
            if not bash:
                continue
            docstr = subprocess.getoutput(f'{command} -h')
        else:
            if not py:
                continue
            docstr = subprocess.getoutput(f'{command} --fullhelp')
        print(f"Building Docs for {command.resolve()}")

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
        docstr = filter('MissingEnvVarWarning',docstr)
        docstr = filter('DeprecationWarning',docstr)
        docs['cmd'].append({'name':command.name, 'doc': docstr})

    if len(docs['cmd']) == 0:
        docs.pop('cmd')
    if idl:
        docs['idl'] = []
        for command in ['spreduce2d.pro','rm_combine_script.pro',
                        'spreduce1d_empca.pro','spcalib_qa.pro',
                        'spspec_target_merge.pro']:
            pf = list(prodir.glob(f"*/{command}"))
            if len(pf) == 0: continue
            print(f'Building Docs for {command}')

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

    with open(docdir/'doc_legacy.rst', 'w') as out:
        out.write(':tocdepth: 2\n\n')
        out.write('.. highlight:: none\n\n')
        out.write(headline('Full Legacy Command Documention') + '\n')
        out.write('Documented below are the legacy CLI primary commands used to run the BOSS Data Reduction Pipeline. '+
                  'However, there are numerous other routines included in this package, which are called by these commands '+
                  'and have their own internal documentation. Additionally, these have been replaced with a updated :doc:`CLI<doc>`\n\n')
        if 'cmd' in docs:
            sec_hdr = 'Full {ctype} Command Usage'
            ctype = []
            if bash:
                ctype.append('Bash')
            if py:
                ctype.append('Python')
            ctype = ' and '.join(ctype)
            sec_hdr = sec_hdr.format(ctype=ctype)
            out = sec(out, sec_hdr, docs['cmd'], 'bin')

        if 'idl' in docs:
            sec_hdr ='IDL Command Usage'
            out = sec(out, sec_hdr, docs['idl'], 'idl')

        out.write('\n.. highlight:: defaults\n\n')

        out.write('\n.. End of document\n')



if __name__ == '__main__' :
    """
    Build BOSS DRP Documention
    """
    build_legacy_docs()
