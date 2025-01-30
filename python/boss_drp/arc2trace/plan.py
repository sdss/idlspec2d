
from pydl.pydlutils.yanny import read_table_yanny
import os
import os.path as ptt


class Plan:
    def __init__(self, mjd, obs, cams = None, planfile=None, outdir = None, vers=None):
        self.mjd = mjd
        self.obs = obs
        self.traceflat = None
        self.tracearc = None
        self.arcs = None
        self.flats = None
        
        self.cams2channel = {'b1':0,'r1':1,'b2':2,'r2':3}
        if obs == 'lco' :
            plan.data_env = 'BOSS_SPECTRO_DATA_S'
            if cams is None : self.cams = ['b2','r2']
            elif isinstance(cams,str) : self.cams = [cams]
        else :
            plan.data_env = 'BOSS_SPECTRO_DATA_N'
            if cams is None : self.cams = ['b1','r1']
            elif isinstance(cams,str) : self.cams = [cams]
        self.channels = [cams2channel[cam] for cam in cams]
    
        
        if vers is None:
            vers = os.environ['RUN2D']
        if outdir is None:
            self.outdir = ptt.join(os.environ['BOSS_SPECTRO_REDUX'],vers, 'trace', f'{mjd}')
        else:
            self.outdir = outdir
            
        self.outhtml = ptt.join(outdir, f'arcs_{mjd}_{obs.lower()}')
        if planfile is None:
            self.planfile = ptt.join(self.outdir,f'spPlanTrace-{mjd}_{obs.upper()}.par')
        else:
            self.planfile = planfile
        if not ptt.exists(self.planfile) :
            raise FileNotFoundError(f'no planfile: {self.planfile}')
        
    def read(self):
        plan = read_table_yanny(self.planfile, 'SPEXP')
        plan.convert_bytestring_to_unicode()
        self.traceflat=plan[plan['flavor'] == 'TRACEFLAT']['name'][0]
        self.tracearc=plan[plan['flavor'] == 'TRACEARC']['name'][0]
        self.arcs=plan[plan['flavor'] == 'arc']
        self.flats=plan[plan['flavor'] == 'flat']

    def flat2traceflat(self, channel):
        flatfile = self.raw2out(self.traceflat, channel)
        flatfile = flatfile.replace('sdR','spTraceFlat')
        return flatfile

    def arc2tracetab(self, arcfile, channel):
        arcfile = self.raw2out(arcfile, channel, exten='.fits.gz')
        arcfile = spTraceTab.replace('sdR','spTraceTab')
        return arcfile

    def raw2out(self, file, channel, exten='.fits.gz'):
        file = file[channel%2].astype(str)
        file = file.replace('.gz','').replace('.fits','').replace('.fit','')
        return file+exten
