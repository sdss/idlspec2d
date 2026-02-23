from boss_drp.utils.getcard import getcard
from boss_drp import MissingEnvVarError
from boss_drp.field import Fieldtype
from boss_drp.utils.splog import splog

try:
    from sdss_access.path import Path
    from sdss_access import Access
    sdss_access = True
except:
    sdss_access = False
from os import getenv
import os.path as ptt
from sys import exit
from pydl.pydlutils.yanny import read_table_yanny, yanny, write_table_yanny

class Sphdrfix:
    def __init__(self, mjd, fps=None, obs='APO', release=None, no_remote=True):
        self.mjd = str(mjd)
        self.sphdrfix_table = None
        if fps is None:
            ft = Fieldtype(mjd=int(mjd))
            fps = ft.fps
        self.fps = fps
        self.obs = obs.lower()
        self.release = release
        if not sdss_access:
            self.release = None
        self.no_remote = no_remote
        self.read_sphdrfix()
        
    def read_sphdrfix(self):
        if self.release is not None:
            path = Path(release=self.release, preserve_envvars=True)
            path_options = {'mjd':self.mjd}
            if path.exists('sdHdrFix', **path_options):
                reportfile = path.full('sdHdrFix', **path_options)
            elif path.exists('sdHdrFix', **path_options, remote=(not self.no_remote)):
                access = Access(release=self.release)
                reportfile = path.full('sdHdrFix', **path_options)
                access.remote()
                access.add('sdHdrFix', **path_ops)
                access.set_stream()
                valid = access.commit()
                if valid is False:
                    return
            else:
                return
        elif not self.fps:
            speclog_dir = getenv('SPECLOG_DIR')
            if speclog_dir is None:
                if splog is None:
                    raise MissingEnvVarError('ERROR: Must set environment variable SPECLOG_DIR')
                else:
                    splog.info('ERROR: Must set environment variable SPECLOG_DIR')
                    exit(1)
            reportfile = ptt.join(speclog_dir, self.mjd, 'sdHdrFix-'+self.mjd+'.par')
        else:
            speclog_dir = getenv('SDHDRFIX_DIR')
            if speclog_dir is None:
                if splog is None:
                    raise MissingEnvVarError('ERROR: Must set environment variable SDHDRFIX_DIR')
                else:
                    splog.info('ERROR: Must set environment variable SDHDRFIX_DIR')
                exit(1)
            reportfile = ptt.join(speclog_dir, self.obs, 'sdHdrfix','sdHdrFix-'+self.mjd+'.par')
        if ptt.exists(reportfile):
            self.sphdrfix_table = read_table_yanny(reportfile, 'OPHDRFIX')
            self.sphdrfix_table.convert_bytestring_to_unicode()

    def fix(self, infile, hdr):
        fileroot = ptt.basename(infile).split('.')[0]
        wfileroot = fileroot.split('-')
        wfileroot = '-'.join([wfileroot[0], '??', wfileroot[-1]])
        if self.sphdrfix_table is not None:
            for row in self.sphdrfix_table:
                if (row['fileroot'] == fileroot) or (row['fileroot'] == wfileroot):
                    hdr[row['keyword']] = row['value']
        if getcard(hdr,'QUALITY') is None:
            hdr['QUALITY'] = 'excellent'
        return hdr
