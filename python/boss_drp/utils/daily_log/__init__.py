from .html import daily_log_html
from .file import daily_log_to_file
from .file import daily_log_js as js
from .email import daily_log_email
from .index import daily_log_index
from .parse_log import (Crash_log, errors, py_err, noerr_cal_b, noerr_cal_r,
                        parse_log, LogCheck, CheckRedux)
from .summary import (_Summary as summary, _summary_html, trace as trace_summary)



def valid_mjd(mjd, mjd_arg, mjdstart, mjdend):
    if mjd_arg is not None:
        return mjd == mjd_arg
    if mjdstart is not None and mjd < mjdstart:
        return False
    if mjdend is not None and mjd > mjdend:
        return False
    return True