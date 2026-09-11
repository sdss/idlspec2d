import re

class Flag:
    def __init__(self,name,color,code,desc):
        self.name=name
        self.color=color
        self.code=code
        self.desc=desc
    def key(self):
        return(f"      <TR><TD><p style='color:{self.color};'>{self.name}</p><TD>{self.desc}</TR>\n")
        
incomplete = Flag('MAGENTA','magenta','#FF00FF','Reduction has yet to be started due to incomplete transfer')
stopped = Flag('RED','red','#FF0000','Stopped Reduction')
NoExp = Flag('MAROON','Maroon','#FFA500','Stopped Reduction for NO GOOD EXPOSURES')
Error_warn = Flag('ORANGE','DarkOrange','#FF8C00','Pipeline ran with errors/warnings')
running = Flag('YELLOW','Gold', '#FFD700','Pipeline is still running')
#NoRedux = Flag('BLUE','blue','#0000FF','No Reductions')
NoObs = Flag('BLUE','blue','#0000FF','No Observations')
NoRedux = Flag('TEAL','MediumTurquoise','#48D1CC','No Science Observations')
#NoObs = Flag('TEAL','Teal','#008080','No Observations')
NoIssues = Flag('GREEN','green','#008000','No issues')
Silent_warn = Flag('ORANGE','#FF8C00','#FF8C00','Pipeline ran with acceptable warnings')


class Complete:
    def __init__(self, custom = None):
        self.line = {}
        self.set(custom)

    def set(self, custom):
        if custom is None:
            self.line = {'spfibermap':'Successful completion of readfibermaps',
                        'spDiag2d':'Successful completion of SPREDUCE2D',
                        'spDiagcomb':'Successful completion of SPCOMBINE',
                        'spDiag1d':'Successful completion of SPREDUCE1D',
                        'spXCSAO':'CPU time to compute RVs',
                        'fieldlist':'Successful completion of fieldlist',
                        'spAll':'Successful completion of fieldmerge',
                        'spSpec_reformat':'Successful completion of spSpec_reformat',
                        'spCalib_QA':'SpectroPhoto QA Complete',
                        'run_spTrace':'Successful completion of boss_arcs_to_trace'}
        else:
            self.line = {'spDiagcomb':'SPSPEC_TARGET_MERGE: Successful completion of spspec_target_merge',
                         'spDiag1d':'Successful completion of SPREDUCE1D',
                         'spXCSAO':'CPU time to compute RVs',
                         'spAll':'Successful completion of fieldmerge',
                         'spSpec_reformat':'Successful completion of spSpec_reformat'}


complete = Complete()


class Crash_log:
    def __init__(self, step, error,msg=None,line=None, flag=Error_warn):
        self.step  = step
        self.error = re.compile('.*'+error+'.*')
        self.msg = msg
        self.line = line
        self.flag = flag
    def check(self,i, line, step):
        if self.step is not None:
            if self.step != step:
                return
        if self.line is not None:
            if i > self.line:
                return
        if self.error.match(line):
            if self.msg is None:
                return(line.replace('\n',''))
            else:
                return(self.msg.format(step=step))
   
   
errors = [Crash_log('spfibermap',' Warning: No matching Field found for DesignID',
                    msg='No matching Field found for DesignID', flag=Error_warn),
          Crash_log('spfibermap',' Warning: No Design Mode found for DesignID',
                    msg='Warning: No Design Mode found for DesignID', flag=Error_warn),
          Crash_log('spfibermap','Warning: SDSS_IDs not found for .* science targets',
                    msg='Warning: SDSS_IDs not found for some science targets', flag=Error_warn),
          Crash_log('spDiag2d','LOCATESKYLINES:.*WARNING: Maximum sky-line shift is.*(DISABLING)'),
          Crash_log('spDiag2d','ABORT: Only            0 sky fibers found',
                    msg='No Sky Fibers Found', flag=stopped),
          Crash_log('spDiag2d','ABORT: No good flats (saturated?)', flag=stopped),
          Crash_log('spDiag2d','SPCALIB: .*: .* paired with no arc', flag=stopped),
          Crash_log('spDiag2d','SUPERFLAT: .*: Creating superflat from .* fibers',
                    flag=stopped, line =1),
          Crash_log('spDiag2d','ABORT: Reject science as too bright: 25-th-percentile =',
                    msg='Reject Bright Science', flag=stopped),
          Crash_log('spDiag2d','SKYSUBTRACT:.*: Discarding .*(fractional) of the sky pixels as bad',
                    msg='Failed Sky Subtraction', line = -1, flag=stopped),
          Crash_log('spDiag2d','FITSPECTRARESOL: .*: Calculating the spectra resolution',
                    msg='Failed FITSPECTRARESOL', line = 1, flag=stopped),
          Crash_log('spDiag2d','EXTRACT_BUNDLE_IMAGE: .*: sigmasize:',
                    msg='Failure Extracting Exposure', line = 1, flag=stopped),
          Crash_log('spDiag2d','FITMEANX: .*:',msg='Failure in Sky Line Identification',
                    line = 1, flag=stopped),
          Crash_log('spDiag2d','ABORT: ERROR in MATCH_TRACE: CHOLDC: choldc failed.', flag=stopped),
          Crash_log('spDiag2d','XCEN is not sorted or not separated by greater than 3 pixels.',
                    msg='Warning: Close or Overlapping Traces', flag=Error_warn),
          Crash_log('spDiag2d','Big wavelength gap',flag=Silent_warn),
          Crash_log('spDiagcomb','RM_SPFLUX_V5:.*: USING XYFIT', flag=stopped,
                    msg='SpectroPhoto Calibration Failure', line = 1),
          Crash_log('spDiagcomb','RM_SPCOMBINE_V5: ABORT: No exposures with SCORE > 0',
                    msg='No Good Exposures', flag=NoExp),
          Crash_log('spDiagcomb','RM_SPFLUX_V5: Rejected  .* of  .* std stars',
                    msg='Failure Combining Exposures', line = 1, flag=stopped),
          Crash_log('spDiagcomb','RM_SPFLUX_V5: ABORT: No good fluxing stars!',
                    flag=Error_warn, msg='ABORT: No good fluxing stars!'),
          Crash_log('spDiagcomb','RM_SPFLUX_V5: WARNING: Already rejected .* of  .* std stars',
                    flag=Error_warn),
          Crash_log('spDiagcomb','RM_SPFLUX_V5: Iteration #',
                    msg='Failure in Fluxing', line = 1, flag=stopped),
          Crash_log('spDiag1d','ZCOMPUTE: .*',msg='Failure in COMPUTECHI2 for ZFIND',
                    line = 1, flag=stopped),
          Crash_log('spDiag1d','ZFIND: .*',msg='Failure in COMPUTECHI2 for ZFIND',
                    line = 1, flag=stopped),
          Crash_log('run_spTrace','Execution halted', msg='Failed run_spTrace', flag=stopped),
          Crash_log('run_spTrace','Killed', msg='Failed run_spTrace', flag=stopped),
          Crash_log('spAll','fieldmerge: EXITING!!', flag=stopped),
          Crash_log('spSpec_reformat', 'read_spAll: ERROR: Missing .*',
                    msg='Failed spSpec_reformat: missing spAll field', flag=stopped)]

py_err = [Crash_log(None,'exception:',
                     msg='Failed {step}', flag=stopped),
         Crash_log(None,'SyntaxError:',
                     msg='Failed {step}', flag=stopped),
         Crash_log('spAll','fieldmerge: No valid spAll entries', flag=stopped),
         Crash_log(None,'FileNotFoundError', msg='Failed {step}', flag=stopped)]

noerr_cal_b = Crash_log('spDiag2d','SPCALIB: b.*: .* paired with arc',
                        msg='No error', flag=NoIssues)
noerr_cal_r = Crash_log('spDiag2d','SPCALIB: r.*: .* paired with arc',
                        msg='No error', flag=NoIssues)
