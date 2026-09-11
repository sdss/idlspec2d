#!/usr/bin/env python3

import faulthandler
from boss_drp import MOUNTAIN 
import os

if MOUNTAIN:
    faulthandler.enable()  # Dumps a traceback on segmentation fault

### These are now imported lazily in the functions that need them
#from boss_drp.sos.SOS import SOS
#from boss_drp.sos import plot
#from boss_drp.sos.db.plot_robodamus import plot_robodamus, robodamus
#from boss_drp.utils.hash import create_hash, check_hash
#from boss_drp.sos import log2html
#from boss_drp.sos.db.loadSN2Value import loadSN2Values
#from boss_drp.sos import parse_runtime
#from boss_drp.sos.build_combined_html import build_combine_html
#from boss_drp.sos.read_sos import read_SOS
#from boss_drp.sos.BOSS_log import build_log
#####

from boss_drp.cli.boss_drp.tools import run_sdR_hdrfix, run_flag_manual_cal
from boss_drp.cli.boss_drp.run import run_boss_arcs_to_trace as _run_boss_arcs_to_trace
from boss_drp.sos.sos_classes import SOS_config
from boss_drp.utils import jdate
from boss_drp.utils.argparse_help import full_help_callback, OrderedGroup
from astropy.time import Time
import numpy as np
import time
import sys
import click
from multiprocessing import Process
import re

def parse_num_list(ctx, param, value):
    """
    Replace this with your real parseNumList logic.
    Return a list of ints or whatever your app expects.
    """
    if value is None:
        return [None]

    m = re.match(r'(\d+)(?:-(\d+))?$', value)
    # ^ (or use .split('-'). anyway you like.)
    if not m:
        raise click.BadParameter("'" + value + "' is not a range of number. Expected forms like '0-5' or '2'.", 
                           ctx=ctx, param=param)
    start = int(m.group(1), 10)
    end = int(m.group(2) or m.group(1), 10)
    return list(range(start, end + 1))

def require_exactly_one(group_name, values):
    selected = [name for name, val in values.items() if val]
    if len(selected) != 1:
        raise click.UsageError(
            f"Exactly one of {group_name} must be selected: {', '.join(values.keys())}"
        )
    return selected[0]


@click.group(
    cls = OrderedGroup,
    invoke_without_command=True,
    context_settings={"help_option_names": ["-h", "--help"],"max_content_width": 150},
)

@click.option("-r","--red", "CCDs", flag_value="red", default=None, help="Red Camera Process")
@click.option("-b","--blue", "CCDs", flag_value="blue", default=None, help="Blue Camera Process")
@click.option("-j","--joint", "CCDs", flag_value="joint", default=None, help="Both Camera Processes")

@click.option("-c","--catchup", "mode", flag_value="catchup", default=None, help="Run Catchup on the night or (MJD)")
@click.option("-t","--redoMode", "mode", flag_value="redoMode", default=None, help="Save outputs of MJD or exposure to sosredo")
@click.option("-d","--test", "mode", flag_value="test", default=None, help="Save outputs and logs to sosredo/dev")
@click.option("--utah", "mode", flag_value="utah", default=None, hidden=True)
@click.option("--systemd", "mode", flag_value="systemd", default=None, hidden=True)

@click.option("--unlock", is_flag=True, default=False, help="Unlock Locked Files")
@click.option("-e", "--exp", callback=parse_num_list, default=None, metavar='EXPID',
               help="exposure id (or range of exp id 500-510) (with or without leading zeros)")
@click.option("-m", "--mjd", multiple=True, type=str, help="MJD", metavar="MJD")
@click.option("--apo", is_flag=True, default=False, hidden=True)
@click.option("--lco", is_flag=True, default=False, hidden=True)
@click.option("--nodb", is_flag=True, default=False, help="skip opsdb load")
@click.option("--no-gz", 'no_gz', is_flag=True, default=False, help="Overrides the requirement for '.gz' compressed files (experimental)")
@click.option("--no-reject",'no_reject', is_flag=True, default=False, help="Overrides the Calibration rejection (use with caution)")
@click.option("--no-flagbad","run_flag_bad", is_flag=True, default=True, help="Update Bad Calibration headers with QUALITY=bad")
@click.option("-f", "--clobber_fibermap", is_flag=True, default=False, help="Clobbers the existing spfibermap files")

@click.option("sdssv_sn2","--no-sdssv-sn2",is_flag=True, default=False, help="Report a second set of SN2 values with updated fit parameters")
@click.option("--no-sn2-15", "sn2_15", is_flag=True,default=True, help="Skip reporting a set of SN2 values with a fiducial mag of 15 for engineering fields")
@click.option("--bright", is_flag=True, default=False, help="Display BOSS_only Bright Time Operation SN2_15 for all fields")

@click.option("-n", "--no-arc2trace", "arc2trace",is_flag=True, default=True, help="Skip Utilizing arc2trace refinements")
@click.option("-o", "--forcea2t", is_flag=True, default=False, help="Force arc2trace for all fields (even if flat exists for field)")
@click.option("--plot", is_flag=True, default=False, help="Produce Science Plots for each exposure")
@click.option("-v", "--verbose", is_flag=True, default=False, help="prints the only (or red if joint) active SOS process to terminal")
@click.option(
    "--fullhelp",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=full_help_callback,
    help="Show full help including all subcommands"
)
@click.pass_context
def cli(ctx, CCDs, mode, unlock, exp, mjd, apo, lco, nodb, no_gz, no_reject,
        clobber_fibermap, sdssv_sn2, sn2_15, bright, arc2trace,
        forcea2t, plot, verbose):
    """SOS process for reducing BOSS data on the Moutain"""
    ctx.ensure_object(dict)
    ctx.obj.update(
        CCDs=CCDs,
        mode=mode,
        unlock=unlock,
        exp=exp,
        mjd=list(mjd) if mjd else None,
        apo=apo,
        lco=lco,
        nodb=nodb,
        no_gz=no_gz,
        no_reject=no_reject,
        clobber_fibermap=clobber_fibermap,
        sdssv_sn2=sdssv_sn2,
        sn2_15=sn2_15,
        bright=bright,
        arc2trace=arc2trace,
        forcea2t=forcea2t,
        plot=plot,
        verbose=verbose,
    )

    # If no subcommand is given, run the legacy behavior.
    if ctx.invoked_subcommand is None:
        run_sos(ctx.obj)


def run_sos(args):
    from boss_drp.sos.SOS import SOS

    require_exactly_one(
        "mode",
        {
            "--catchup": args["mode"] == "catchup",
            "--redoMode": args["mode"] == "redoMode",
            "--test": args["mode"] == "test",
            "--utah": args["mode"] == "utah",
            "--systemd": args["mode"] == "systemd",
        },
    )

    if args["apo"]:
        os.environ["OBSERVATORY"] = "APO"
    if args["lco"]:
        os.environ["OBSERVATORY"] = "LCO"

    obs = os.getenv("OBSERVATORY")
    if obs and obs.upper() == "APO":
        blue = "b1"
        red = "r1"
    else:
        blue = "b2"
        red = "r2"

    ccd_mode = args["CCDs"]

    require_exactly_one(
        "CCDs",
        {
            "--red": args["CCDs"] == "red",
            "--blue": args["CCDs"] == "blue",
            "--joint": args["CCDs"] == "joint",
        },
    )
    if ccd_mode == "red":
        args["CCDs"] = [red]
    elif ccd_mode == "blue":
        args["CCDs"] = [blue]
    elif ccd_mode == "joint":
        args["CCDs"] = [blue, red]


    if not args["arc2trace"]:
        args["forcea2t"] = False

    mjds = args["mjd"] or [None]

    for _mjd in mjds:
        procs = {}
        for i, CCD in enumerate(args["CCDs"]):
            if args["CCDs"] in ([red], [blue]):
                pause = False
                termverbose = args["verbose"]
            elif CCD == red:
                pause = True
                termverbose = True if args["verbose"] else False
            else:
                pause = True
                termverbose = False

            kwrds = dict(
                exp=args["exp"],
                mjd=_mjd,
                catchup=(args["mode"] == "catchup"),
                redoMode=(args["mode"] == "redoMode"),
                systemd=(args["mode"] == "systemd"),
                nodb=args["nodb"],
                no_reject=args["no_reject"],
                clobber_fibermap=args["clobber_fibermap"],
                sdssv_sn2=args["sdssv_sn2"],
                sn2_15=args["sn2_15"],
                arc2trace=args["arc2trace"],
                forcea2t=args["forcea2t"],
                bright_sn2=args["bright"],
                pause=pause,
                test=(args["mode"] == "test"),
                unlock=args["unlock"],
                plot=args["plot"],
                utah=(args["mode"] == "utah"),
                termverbose=termverbose,
            )

            if len(args["CCDs"]) == 1:
                SOS(CCD, **kwrds)
            else:
                procs[i] = Process(target=SOS, args=(CCD,), kwargs=kwrds)
                procs[i].start()

        if len(args["CCDs"]) > 1:
            try:
                for i in range(len(args["CCDs"])):
                    procs[i].join()
            except KeyboardInterrupt:
                print("Main process interrupted. Terminating child processes...")
                for i in range(len(args["CCDs"])):
                    procs[i].terminate()
                sys.exit(1)

@cli.group(name='Tools', context_settings={"help_option_names": ["-h", "--help"],"max_content_width": 150})
@click.option(
    "--fullhelp",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=full_help_callback,
    help="Show full help including all subcommands"
)
def tools():
    """Tools to use with SOS"""
    pass

@tools.command(name='plot', context_settings={"help_option_names": ["-h", "--help"],"max_content_width": 150})
@click.option('-m','--mjd', help='MJD of reduction', type=str, required=True, metavar='MJD')
@click.option('-e','--expid', help='Exposure ID to plot', type=str, required=True, metavar='EXPID')
@click.option('-o','--observatory', "obs", envvar="OBSERVATORY", metavar='OBS',
              help='Observatory (default: $OBSERVATORY)')
@click.option('--ccd', multiple=True, default=None, help=f'CCDs to plot; defaults to both CCDs')
@click.option('--redo', is_flag=True, default=False,  help='If set use sosredo rather then sos reductions')
@click.option('--dev', is_flag=True, default=False,  help='If set use sosredo/dev rather then sos reductions')
@click.option('--mask_end', is_flag=True, default=False,  help='Mask end of the spectra during plotting')

@click.option('--ToOs', is_flag=True, default=False, help='Plot only ToO fibers')
@click.option('--assigned', is_flag=True, default=False, help='Plot only fibers assigned to targets (includes ToOs)')
@click.option('--science', is_flag=True, default=False,  help='Plot all fibers with science targets (includes ToOs and assigned)')
@click.option('--pdf', is_flag=True, default=False,  help='Plot all fibers into a single multi-panel PDF instead of individual PNGs')
@click.option('--single_ccd', is_flag=True, default=False, help='Only plot the specified CCD even if both are available')
def run_plot(mjd, expid, obs, ccd, redo, dev, mask_end, toos, assigned, science, pdf, single_ccd):
    """Plot the Science frame for SOS"""
    from boss_drp.sos import plot
    if not ccd:
        ccd = ['b2','r2'] if obs == 'LCO' else ['b1','r1']
    if single_ccd: 
        ccd = [ccd[0]]
    for _ccd in ccd:
        plot(mjd, expid, obs, _ccd, redo=redo, dev=dev, mask_end=mask_end, 
             ToOs=toos, assigned=assigned, science=science, pdf=pdf, single_ccd=single_ccd)

@tools.command(name='plot_night', context_settings={"help_option_names": ["-h", "--help"],"max_content_width": 150})
@click.option('-m','--mjd', help='MJD of reduction', type=str, required=True, metavar='MJD')
@click.option('-o','--observatory', "obs", envvar="OBSERVATORY", metavar='OBS',
              help='Observatory (default: $OBSERVATORY)')
@click.option('--ccd', multiple=True, default=None, help=f'CCDs to plot; defaults to both CCDs')
@click.option('--redo', is_flag=True, default=False,  help='If set use sosredo rather then sos reductions')
@click.option('--dev', is_flag=True, default=False,  help='If set use sosredo/dev rather then sos reductions')
@click.option('--mask_end', is_flag=True, default=False,  help='Mask end of the spectra during plotting')

@click.option('--ToOs', is_flag=True, default=False, help='Plot only ToO fibers')
@click.option('--assigned', is_flag=True, default=False, help='Plot only fibers assigned to targets (includes ToOs)')
@click.option('--science', is_flag=True, default=False,  help='Plot all fibers with science targets (includes ToOs and assigned)')
@click.option('--pdf', is_flag=True, default=False,  help='Plot all fibers into a single multi-panel PDF instead of individual PNGs')
@click.option('--single_ccd', is_flag=True, default=False, help='Only plot the specified CCD even if both are available')
def run_plot(mjd, obs, ccd, redo, dev, mask_end, toos, assigned, science, pdf, single_ccd):
    """Plot the Science frame for SOS"""
    from boss_drp.sos import plot
    if not ccd:
        ccd = ['b2','r2'] if obs == 'LCO' else ['b1','r1']
    if not single_ccd: 
        ccd = [ccd[0]]
    for _ccd in ccd:
        plot(mjd, None, obs, _ccd, redo=redo, dev=dev, mask_end=mask_end, 
             ToOs=toos, assigned=assigned, science=science, pdf=pdf, single_ccd=single_ccd)

@tools.command(name='robodamus', context_settings={"help_option_names": ["-h", "--help"],"max_content_width": 150})
@click.option('--mjd','-m', help='SJD of reduction', type=float, default=Time.now().mjd + 0.3, metavar='MJD')
def run_robodamus(mjd):
    """Plot the Robodamus predictions vs the SOS SN2 Measurements"""
    from boss_drp.sos.db.plot_robodamus import plot_robodamus, robodamus
    if not MOUNTAIN:
        raise click.BadOptionUsage('robodamus is not configured to run off of the mountains...')

    if mjd < 0:
        mjd = jdate.sjd()+mjd
    mjd = int(np.floor(mjd))
    output = robodamus(sjd = mjd)
    dest = '/data/boss/sos/tests/robodamus' 
    if not os.path.exists(dest):
        os.makedirs(dest, exist_ok=True)
    plot_robodamus(output, os.path.join(dest,f"robodamus_sn_{int(np.floor(mjd))}.png"))

@tools.command(name='hash', context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 150}) 
@click.option("--mjd", type=int, help="MJD to process",required=True, metavar='MJD')
@click.option("--redo",'-t', is_flag=True, default=False, help='use sosredo directory (same as SOS -t or --redoMode option)')
@click.option("--test",'-d', is_flag=True, default=False, help='use sosredo/dev directory (same as SOS -d or --test option)')
@click.option("--utah",'-u', is_flag=True, default=False, help='use utah test directory (same as SOS --utah option)')
@click.option("--create", is_flag=True, default=False, help="Create the Hash file")
@click.option("--check", is_flag=True, default=False, help="Check the SOS Hash")
@click.option("--transfer",is_flag=True, default=False, help="Check the SOS Hash of a Utah transfer")
@click.option("--lco", is_flag=True, default=False, help="Build/check for lco at Utah")
@click.option("--dummy", is_flag=True, default=False, help="Create a dummy file to prevent an empty hash file")
def run_hash(mjd, redo, test, utah, create, check, transfer, lco, dummy):
    """Create or check the SOS file hash"""
    from boss_drp.utils.hash import create_hash, check_hash
    if transfer:
        sosdir = os.getenv('BOSS_SOS_S') if lco else os.getenv('BOSS_SOS_N')
        if not sosdir:
            raise click.ClickException("Missing BOSS_SOS_S/BOSS_SOS_N for --transfer")

    else:
        if utah:
            os.environ['OBSERVATORY'] = 'APO' if not lco else 'LCO'
        SOS_config.setup(mjd=mjd, redo=redo, test=test, utah=utah)
        sosdir = SOS_config.sosdir

    sosdir = os.path.join(sosdir, f'{mjd}')
    if create:
        i = 0
        while create_hash(sosdir, dummy=dummy):
            i += 1
            if i >= 10:
                print("Error: Could not unlock sha1sum file")
                break
            time.sleep(3)
            print("\nsha1sum is locked")

    if check:
        check_hash(sosdir, verbose = True)

@tools.command(name='log2html', context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 150}) 
@click.option('--mjd', help='MJD of reduction', type=str, required=True, metavar='MJD')
@click.option('--sosdir', type=str, required=True, metavar='SOSDIR',
              help='Path to SOS save directory (ie. folder that contains logfile-?????.fits)')
@click.option('-l','--logfile', help='Name of logfile (default: logfile-{mjd}.fits)', metavar='LOGFILE')
@click.option('-f','--htmlfile', help='Name of output htmlfile (default: logfile-{mjd}.html)', metavar='HTMLFILE')
@click.option('--obs','-o', help=f"Observatory of observations (default: {os.getenv('OBSERVATORY','apo')})", metavar='OBS')
@click.option('--copydir', '-c', help='Where to save the htmls', default=None, metavar='SAVEDIR')
@click.option('--fps', is_flag=True, default=False,  help='build for FPS reductions')
@click.option('--sdssv_sn2', is_flag=True, default=False,  help='Include SDSSV SN2 V2')
@click.option('--sn2_15', is_flag=True, default=False,  help='Include Mag 15 SN2')
@click.option('--bright', is_flag=True, default=False,  help='Include Mag 15 SN2 for all exposures')
def run_log2html(mjd, sosdir, logfile, htmlfile, obs, copydir, fps, sdssv_sn2, sn2_15, bright):
    """Create the HTML Logging Page for SOS"""
    from boss_drp.sos import log2html

    log2html(mjd, sosdir, logfile=logfile, htmlfile=htmlfile,
             obs = obs, fps=fps, sdssv_sn2=sdssv_sn2,
             sn2_15 = sn2_15, bright = bright, copydir = copydir)


@tools.command(name='loadsn2', context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 150}) 
@click.option('--fits',help='The fits file is the science frame output from sos-reduce', required=True, metavar='FITSFILE')
#@click.option('--confSum',help='confSummary-file', required=True, metavar='CONFSUMMARY')
@click.option('-v','--verbose', is_flag=True, default=False, help='verbose')
@click.option('-u','--update',  is_flag=True, default=False, 
              help='update (An error will occur if the exposure has already been processed, unless set)')
@click.option('--sdssv_sn2',    is_flag=True, default=False, help='Load sdssv_sn2')
def run_loadsn2(fits, verbose, update, sdssv_sn2):
    """Load SOS SN2 values into OpsDB"""
    from boss_drp.sos.db.loadSN2Value import loadSN2Values
    loadSN2Values(fits, verbose=verbose, update=update, sdssv_sn2=sdssv_sn2)


@tools.command(name='parse_runtime', context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 150}) 
@click.argument('LogFile')
@click.option('-a','--all', is_flag=True, default=False, help='Combine all daily logs of this format')
@click.option('-s','--stamp',is_flag=True, default=False, help='Add Date Stamp to output file')
def run_parse_runtime(logfile,all,stamp):
    """Process log file to calculate elapsed times for SOS"""
    from boss_drp.sos import parse_runtime
    parse_runtime(logfile,all=all,stamp=stamp)


@tools.command(name='htmlIndex', context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 150}) 
@click.option('-s','--sosdir', help="Base SOS output directory", required=True, metavar='SOSDIR')
@click.option('-f','--force', is_flag=True, default=False, help='Force Update of Index page')
def run_htmlIndex(sosdir,force):
    """Build sos/combined/index.html"""
    from boss_drp.sos.build_combined_html import build_combine_html
    build_combine_html(sosdir,force=force)


@tools.command(name='FiberQA', context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 150}) 
@click.option('-s','--sosdir', help="Base SOS output directory",required=True, metavar='SOSDIR')
@click.option('-m','--mjd', help='MJD', required=True, metavar='MJD')
@click.option('-e','--exp', help='Exposure Name', metavar='EXPID')
@click.option('-n','--nocopy', is_flag=True, default=False, help='Prevent copy to combined Directory')
@click.option('--no-hash','no_hash', is_flag=True, default=False, help='Skip updating the file hash')
def run_FiberQA(sosdir,mjd,exp,nocopy,no_hash):
    """Create Fiber info Summary for SOS"""
    from boss_drp.sos.read_sos import read_SOS
    read_SOS(sosdir, mjd, exp, nocopy=nocopy, update_hash=(not no_hash))



def add_observatory_option(f):
    if not MOUNTAIN:
        f = click.option(
            "-o", "--observatory", "--obs",
            default=None,
            type=click.Choice(['apo', 'lco'], case_sensitive=False),
            help="Manually set observatory",
        )(f)
    return f


def log_options(f):
    f = click.option("-m", "--mjd", type=str, default=None, help="MJD")(f)
    f = click.option("-y", "--yesterday",is_flag=True,default=False, help="current mjd-1")(f)
    f = add_observatory_option(f)
    f = click.option("-l", "--long","long_",is_flag=True, default=False, help="Long/detailed version of log")(f)
    f = click.option("--new_ref",is_flag=True,default=False,
                     help="Calculate new reference values in fratio and w_shift and show in place of fratio and w_shift")(f)
    f = click.option("-c", "--hide_hart", "--hide_hartmann","hide_hartmann", is_flag=True,default=False,
                             help="Hide cleaned version of Hartmann Logs as a table")(f)
    f = click.option("-r", "--hart_raw",is_flag=True, default=False,
                             help="Print raw form (instead of table form) of Hartmann Logs")(f)
    f = click.option("-e", "--hide_error", is_flag=True, default=False, help="Hide SOS Error and Workings")(f)
    f = click.option("-s", "--hide_summary",is_flag=True, default=False,help="Hide data summary table")(f)
    f = click.option("-t", "--show_ToOs",is_flag=True, default=False,help="Show ToOs")(f)
    f = click.option('--email', hidden=True)(f)
    return f

@tools.command(name='Log', context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 150}) 
@log_options
def run_log(mjd, yesterday, long_, new_ref, hide_hartmann, hart_raw, hide_error, hide_summary,
            show_toos, observatory=None, email=None):
    """Build BOSS Exposure Log"""
    from boss_drp.sos.BOSS_log import build_log
    if MOUNTAIN:
        observatory = None

    try:
        obs = os.getenv('OBSERVATORY').lower()
        if observatory is not None:
            obs = observatory
    except:
        if observatory is None:
            obs = input('Enter Observatory {apo,lco} ')
            obs = obs.lower()
        else:
            obs = observatory.lower()

    if obs not in ['apo', 'lco']:
        raise click.ClickException(f'Error: Invalid observatory ({obs})')

    if mjd is None:
        mjd = jdate.obs(obs).astype(str)
        if yesterday is True: mjd = str(int(mjd)-1)
    if MOUNTAIN:
        datadir = '/data/spectro/'
        sos_dir = '/data/boss/sos/'
    else:
        hide_hartmann = True
        hart_raw = False
        if obs == 'apo':
            datadir =  os.getenv('BOSS_SPECTRO_DATA_N')
            sos_dir =  os.getenv('BOSS_SOS_N')
        else:
            datadir =  os.getenv('BOSS_SPECTRO_DATA_S')
            sos_dir =  os.getenv('BOSS_SOS_S')
    if hart_raw is True:
        hide_hartmann = False
    build_log(mjd, obs, Datadir=datadir, long_log = long_, new_ref = new_ref,
              hart=not hide_hartmann, hart_table = not hart_raw, hide_error=hide_error,
              hide_summary=hide_summary, too=show_toos, sos_dir = sos_dir, email=email)


tools.add_command(run_sdR_hdrfix)
tools.add_command(run_flag_manual_cal)


mountain_opt = click.Option(
    ["--mountain"],
    is_flag=True,
    default=False,
    help="Run arc2trace refinements in mountain mode",
)
redo_opt = click.Option(["--redo",'-t'], is_flag=True, default=None, help='use sosredo directory (same as SOS -t or --redoMode option)')
test_opt = click.Option(["--test",'-d'], is_flag=True, default=None, help='use sosredo/dev directory (same as SOS -d or --test option)')
utah_opt = click.Option(["--utah",'-u'], is_flag=True, default=None, help='use utah test directory (same as SOS --utah option)')

def _boss_arcs_to_trace_wrapper(**kwargs):
    mountain = kwargs.pop("mountain", False)
    redo = kwargs.pop("redo", None)
    test = kwargs.pop("test", None)
    utah = kwargs.pop("utah", None)
    if mountain:
        kwargs['threads'] = 0
        kwargs['obs'] = 'apo' if os.getenv('OBSERVATORY', 'apo').lower() == 'apo' else 'lco'
        kwargs['vers'] = 'sos'
    if redo or test or utah:
        SOS_config.setup(redo=redo, test=test, utah=utah)
        kwargs['sosdir'] = SOS_config.sosdir
        os.environ['BOSS_SPECTRO_REDUX'] = os.path.join(SOS_config.sosdir,f'{kwargs["mjd"]}')

    return _run_boss_arcs_to_trace.callback(**kwargs)

run_boss_arcs_to_trace_ext = click.Command(
    name=_run_boss_arcs_to_trace.name,
    params=list(_run_boss_arcs_to_trace.params) + [
        mountain_opt, redo_opt, test_opt, utah_opt,
    ],
    callback=_boss_arcs_to_trace_wrapper,
    help=_run_boss_arcs_to_trace.help,
    context_settings=_run_boss_arcs_to_trace.context_settings,
)

tools.add_command(run_boss_arcs_to_trace_ext)


if __name__ == "__main__":
    cli()
