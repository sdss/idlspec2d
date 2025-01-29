from boss_drp.sos.arc_to_trace_soshtml import soshtml
from boss_drp.sos.report_err import report
from boss_drp.utils.splog import splog, splog_name
from pyvista import boss

import sys
import matplotlib
import os
import os.path as ptt
import time

class CaptureOutput:
    def __init__(self, logger):
        self.output_list = []  # List to store captured output
        self._original_stdout = None  # To store the original sys.stdout
    def __enter__(self):
        self._original_stdout = sys.stdout  # Save the original sys.stdout
        sys.stdout = self  # Redirect sys.stdout to this instance
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stdout = self._original_stdout  # Restore the original sys.stdout

    def write(self, message):
        self.output_list.append(message)  # Append the message to the list
        self._original_stdout.write(message)  # Print the message to the console

    def flush(self):
        pass  # Implement if flush functionality is needed (e.g., for print calls)


def boss_arcs_to_traces(mjd = None, outdir = None, obs = 'lco', vers = 'master',
                        threads = 8, nskip = 40, cams = None, fitsname = None,
                        designMode = 'uknown', sosdir = None, clobber = False,
                        capture =CaptureOutput, logger = None):
    matplotlib.use('Agg')
    if vers.lower() == 'sos':
        vers = ''
        sos = True
        if sosdir is not None:
            os.environ["BOSS_SPECTRO_REDUX"]  = ptt.join(f'{sosdir}',f'{mjd}')
    else:
        sos = False
    with capture(logger) as captured_output:
        try:
            import logging
            _log = logging.getLogger("astropy")
            _log.setLevel(logging.CRITICAL)
            mjd = int(mjd)
            boss.arc_transform(mjd, obs=obs, clobber=clobber, threads=threads,
                               outdir=outdir, vers=vers, cams=cams)
        except Exception as e:
            import traceback
            print(traceback.format_exc())
            print(type(e).__name__, ":", e)
            if not sos:
                exit()

        if sos:
            fitsname_base = ptt.splitext(ptt.splitext(ptt.basename(fitsname))[0])[0]
            error_message = f"error with {fitsname_base}".lower()
            has_error = any(error_message in line.lower() for line in captured_output.output_list)
            if has_error:
                message = f"BOSS_ARCS_TO_TRACES: {fitsname}: ABORT: Failure in Fitting Arc to Traces (Please retake calibrations for this field)"
                print(message)
                report(fitsname, cams, obs, mjd, message, DESIGNMODE=designMode)
            soshtml(mjd, obs, sosdir)

        else:
            print(f'Successful completion of boss_arcs_to_trace for mjd {mjd} obs {obs} at {time.ctime()}')
