#!/usr/bin/env python3
from boss_drp.prep.readfibermaps.readfibermaps import readfibermaps as pipe_readfibermaps
from boss_drp.spec1d.run_PyXCSAO import run_PyXCSAO
from boss_drp.sos.arc2tracelogger import Logger
from boss_drp.prep.boss_arcs_to_traces import boss_arcs_to_traces
from boss_drp.utils.hash import create_hash
from boss_drp.post.fieldlist import fieldlist
from boss_drp.post.fieldmerge import fieldmerge
from boss_drp.post.build_target_summary import build_target_summary
from boss_drp.post.spSpec_reformat import spSpec_reformat
from boss_drp.post.update_flags import update_Targeting_flags
from boss_drp.post import plot_QA
from boss_drp.post.spcalib_qa import spcalib_qa

from boss_drp.utils.argparse_help import AttrDict, multi_str2bool, multi_str2none, full_help_callback, OrderedGroup

from sdss_access import Access
from pyvista import boss
import astropy.time

import os
import platform
import click
import sys
from datetime import date


@click.group('run', cls = OrderedGroup, context_settings={"help_option_names": ['-h','--help'], "max_content_width": 150}) 
@click.option(
    "--fullhelp",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=full_help_callback,
    help="Show full help including all subcommands"
)
def run():
    """Miscellaneous BOSS DRP Run Steps"""
    pass


SHOW_NO_DB = ("sdss5" in platform.node()) or (os.getenv("IDLSPEC2D_SOS") is not None)


@run.command(name='readfibermap', context_settings={"help_option_names": ["-h", "--help"]},
    help="Produces spfibermap file corresponding to a spplan2d (or single confSummary file for SOS).",
)
@click.option("-p","--spplan2d",help="spplan2d file for idlspec2d run")
@click.option("--topdir", help=("Alternative output directory (defaults to location of spplan2d file or "
                                "/data/boss/sos/{mjd} for SOS)"))
@click.option("-c","--clobber",is_flag=True, help="Overwrites previous spfibermap file")
@click.option("-n","--no-db", 'no_db', is_flag=True, hidden=not SHOW_NO_DB,
              help="Bypasses SDSSDB access and utilizes MOS target files from SDSS-V DR")
@click.option("--fast", is_flag=True,
              help="When using --no-db, streamlines process and only gets parallax from MOS target files")
@click.option("--datamodel",
              help="Supply a datamodel file (defaults to $IDLSPEC2D/datamodel/spfibermap_dm.par or $IDLSPEC2D/datamodel/spfibermap_sos_dm.par for SOS)")
@click.option("-s","--SOS", "sos", is_flag=True, help="Produces spfibermap for SOS")
@click.option("--release", default="sdsswork", show_default=True,
              help=("sdss_access data release (defaults to sdsswork), required if you do not have proprietary access"))
@click.option("--remote",is_flag=True, help="Allow for remote access to data using sdss-access")
@click.option("--v_targ", "V_TARG", default="*", show_default=True, help="SDSS-V MOS Targeting Product Version (for no Database access use)")
@click.option("--confSummary","confSummary", help="confSummary file for SOS (required with --SOS)")
@click.option("--ccd",type=click.Choice(["b2", "r2", "b1", "r1"], case_sensitive=True), help="CCD for SOS")
@click.option("--mjd", type=str, help="MJD of observation")
@click.option("--log", is_flag=True, help="Creates log file in topdir")
@click.option("--lco", is_flag=True, hidden=True)
@click.pass_context
def run_readfibermap(ctx, spplan2d, topdir, clobber, no_db,fast,datamodel,sos,release, 
                 remote,V_TARG, confSummary, ccd, mjd,log,lco):


    if lco:
        os.environ["OBSERVATORY"] = "LCO"

    if release != "sdsswork":
        if release not in Access().get_available_releases():
            raise click.ClickException(f"ERROR: {release} is not a valid release")
    elif remote:
        try:
            Access().remote()
        except Exception:
            raise click.ClickException(
                "ERROR: No netrc file found. see "
                "https://sdss-access.readthedocs.io/en/latest/auth.html#auth"
            )

    if sos:
        if confSummary is None or ccd is None or mjd is None:
            raise click.UsageError("ERROR: --confSummary, --ccd, and --mjd are required with the --SOS option")

        SOS_opts = {"confSummary": confSummary, "ccd": ccd, "mjd": mjd, "log": log, "log_dir": topdir}
    else:
        if spplan2d is None:
            raise click.UsageError("ERROR: --spplan2d is required without the --SOS option")
        SOS_opts = None

    pipe_readfibermaps(spplan2d=spplan2d, topdir=topdir, clobber=clobber, SOS=sos, no_db=no_db, fast=fast,
                  datamodel=datamodel, SOS_opts=SOS_opts, release=release, remote=remote, V_TARG=V_TARG)

@run.command(name='boss_arcs_to_traces', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.option("--mjd", required=True, type=int, help="MJD to process")
@click.option("--outdir", type=str, help="output directory")
@click.option("--obs", type=click.Choice(["lco", "apo"], case_sensitive=False),
              default="lco", show_default=True, help="observatory")
@click.option("--vers", type=str, default="master", show_default=True,
              help="BOSS_SPECTRO_REDUX version")
@click.option("--threads", type=int, default=8, show_default=True, help="number of threads")
@click.option("--cams", type=str, default=None, help="Supply the camera for operation with SOS files")
@click.option("--fitsname", type=str, default=None, help="Supply the FitsName for SOS error reporting")
@click.option("--sosdir", type=str, default=None, help="Base SOS output directory")
@click.option("--clobber", is_flag=True, help="clobber?")
@click.option("--no-hash",'no_hash', is_flag=True, help="Skip updating the file hash")
def run_boss_arcs_to_trace(mjd, outdir, obs, vers, threads, cams, fitsname, sosdir, clobber, no_hash):
    """Routine to transfer trace locations from an initial arc/flat pair to subsequent arc frames for use with the science frames"""
    if vers.lower() == "sos":
        if not fitsname or not sosdir:
            raise click.UsageError("ERROR: --fitsname and --sosdir are required when --vers is sos")

        logfile = os.path.splitext(os.path.splitext(os.path.basename(fitsname))[0])[0] + ".log"
        logfile = logfile.replace("sdR", "spTraceTab")
        logfile = os.path.join(f"{sosdir}", f"{mjd}", "trace", f"{mjd}", logfile)

        log = Logger(logfile)
        log.start(cmd=sys.argv)

        try:
            boss_arcs_to_traces(mjd=mjd, outdir=outdir, obs=obs, vers=vers, threads=threads,
                                cams=cams, fitsname=fitsname, sosdir=sosdir, clobber=clobber)
        finally:
            log.stop()

        if not no_hash:
            test = create_hash(os.path.join(sosdir, str(mjd)))
            if test:
                print("\nsha1sum is locked")
        return

    boss_arcs_to_traces(mjd=mjd, outdir=outdir, obs=obs, vers=vers, threads=threads,
                        cams=cams, fitsname=fitsname, sosdir=sosdir, clobber=clobber)



@run.command(name='run_PyXCSAO', short_help = 'Runs pyXCSAO to determine RVs',
            context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.argument("fitsfile", type=click.Path(exists=True, dir_okay=False, readable=True))
@click.option("-r", "--run1d", envvar='RUN1D', show_default="env: RUN1D", help="run1d name")
@click.option("--epoch", is_flag=True, help="Run for epoch coadds")
@click.option("--custom", type=str, help="Name of custom coadd")
def cmd_run_PyXCSAO(fitsfile, run1d, epoch, custom):
    """
    Run PyXCSAO for a full spField FITS file using the phoenix_full1 template grid.
    The input file can be either a normal FITS or gzipped FITS file.
    """
    run_PyXCSAO(fitsfile, run1d=run1d, epoch=epoch, custom=custom)


@run.command(name='fieldlist', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.option("-c","--create",is_flag=True, help="Create Fieldlist")
@click.option("--topdir",type=str,envvar="BOSS_SPECTRO_REDUX",show_default="env: BOSS_SPECTRO_REDUX",
               help="Optional override value for $BOSS_SPECTRO_REDUX")
@click.option("--run1d", multiple=True, type=str, show_default="env: RUN1D",
              default=lambda: (os.getenv("RUN1D"),) if os.getenv("RUN1D") else (),
              help="Optional override value for $RUN1D")
@click.option("--run2d", multiple=True, type=str, show_default="env: RUN2D",
              default=lambda: (os.getenv("RUN2D"),) if os.getenv("RUN2D") else (),
              help="Optional override value for $RUN2D")
@click.option("--outdir", type=str, default=None,
              help="Optional output directory (defaults to topdir/$RUN2D)")
@click.option("--skipcart", multiple=True,type=str, default=None,help="List of cartridges to skip")
@click.option("--epoch",is_flag=True, help="Produce FieldList for epoch coadds")
@click.option("--abandoned",is_flag=True, help="Produce FieldList for epoch coadds (including abondoned)")
@click.option("--started",is_flag=True, help="Produce FieldList for epoch coadds (including started)")
@click.option("--basehtml",type=str, help="HTML path for figure (defaults relative to topdir)")
@click.option("--logfile", type=str, default=None, help="Manually set logfile (including path)")
@click.option("--debug", is_flag=True, help="Print full python errors instead of simplified logger messages")
@click.option("--noplot", is_flag=True, help="Skip updating the sky plots")
@click.pass_context
def run_fieldlists(ctx, **kwrds):
    """Build/load BOSS Fieldlist"""
    args = AttrDict(ctx.params)
    args.run1d = list(args.run1d) if args.run1d else None
    args.run2d = list(args.run2d) if args.run2d else None
    args.skipcart = list(args.skipcart) if args.skipcart else None
    _ = fieldlist(**args)


@run.command(name='fieldmerge', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.option("--run2d", type=str, envvar="RUN2D", show_default="env: RUN2D",
              help="Optional override value for the environment variable $RUN2D")
@click.option("--indir", type=str, envvar="BOSS_SPECTRO_REDUX", show_default="env: BOSS_SPECTRO_REDUX",
              help="Optional override value for the environment variable $BOSS_SPECTRO_REDUX")
@click.option("--skip_line", is_flag=True, help="Skip the generation of spAllLine.fits")
@click.option("--include_bad", is_flag=True, help="Include bad fields")
@click.option("--legacy", is_flag=True,
              help="Include columns used by SDSS-IV and depreciated in SDSS-V")
@click.option("--skip_specprimary", is_flag=True,
              help="Skip creation of specprimary and associated columns")
@click.option("--update_specprimary", is_flag=True,
              help="Keep existing specprimary and associated columns and only update new row (and their secondaries)")
@click.option("--lite", is_flag=True, help="Produce lite version of spAll file")
@click.option("--include_XCSAO", "XCSAO", is_flag=True, help="Include XCSAO columns")
@click.option("-f", "--field", type=str, default=None, help="Run for a single Field")
@click.option("-m", "--mjd", type=str, default=None, help="Run for a single MJD")
@click.option("--clobber", "--clobber_fmjd", "clobber", is_flag=True, help="Clobber all spAll-field-mjd files")
@click.option("--clobber_mjd", "clobber_mjd", is_flag=True, help="Clobber all spAll-MJD files")

@click.option("--bkup", is_flag=True, help="Backup existing spAll files")
@click.option("--verbose", is_flag=True, help="Log columns not saved")
@click.option("--logfile", type=str, help="Manually set logfile")
@click.option("--epoch", is_flag=True, help="Produce spAll for epoch coadds")
@click.option("--dev", is_flag=True, hidden=True)
@click.option("--programs", multiple=True,
              help="List of programs to include. Repeat the option for multiple values.")
@click.option("--datamodel", type=str,
              help="Supply a spAll datamodel file (defaults to $IDLSPEC2D/datamodel/spall_dm.par)")
@click.option("--line_datamodel", type=str,
              help="Supply a spline datamodel file (defaults to $IDLSPEC2D/datamodel/spzline_dm.par)")
@click.option("--outroot", type=str,
              help="Path and root of filename for output (defaults to spectra/full or summary)")
@click.option("--allsky", is_flag=True, help="Build spAll for Allsky Custom Coadd")
@click.option("--custom", type=str, help="Name of Custom Coadd")
@click.option("--run1d", type=str,envvar="RUN1D",  show_default="env: RUN1D",
              help="Optional override value for the environment variable $RUN1D (only for custom allsky coadds)")
@click.option("--ndays", type=int, default=None,  help="Limit update to last ndays")
@click.option("--freeze_output", is_flag=True, help="Freeze MJD limited parquet files")
@click.option("--update_target_flags", is_flag=True,
              help="Use the spTargeting file to update the summary file to the latest Targeting Flags")
@click.option("--mjdstart", type=int, default=None, help="Limit update to MJD on/after")
@click.option("--mjdend", type=int, default=None, help="Limit update to MJD on/before")
@click.option("--MJD_dir", type=str, default=None,
              help="Location to save the MJD level temporary files (defaults to BOSS_SPECTRO_SCRATCH)")
@click.option("--to_fits", is_flag=True, help="Dump Parquet to fits format")
@click.option("--force", "--force_rebuild", "force_rebuild", is_flag=True,
              help="Rebuild Summary even if nothing changed")
@click.option("--keep_active", is_flag=True, help='Run "touch" on all intermediate files to keep them active')
@click.pass_context
def run_fieldmerge(ctx, **kwrds):
    """Build BOSS spAll Summary Files"""
    args = AttrDict(ctx.params)
    if args.mjdstart is not None:
        todaymjd = int(float(astropy.time.Time( str(date.today())).jd)-2400000.5)
        args.mjdstart = todaymjd-args.mjdstart

    if args.field is not None:
        fieldmerge(**args)
    else:
        args['clobber'] = args.get('clobber_mjd', False)
        build_target_summary(**args)

@run.command(name='spSpec_reformat', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.option("-f", "--field", type=str, required=True, help="Run for a single Field")
@click.option("-m", "--mjd", type=str, required=True, help="Run for a single MJD")
@click.option("--topdir", type=str, envvar='BOSS_SPECTRO_REDUX',
              show_default="env: BOSS_SPECTRO_REDUX",
              help="Optional override value for the environment variable $BOSS_SPECTRO_REDUX")
@click.option("--run2d", type=str, envvar="RUN2D", show_default="env: RUN2D",
              help="Optional override value for the environment variable $RUN2D")
@click.option("--run1d", type=str, envvar="RUN1D", show_default="env: RUN1D",
              help="Optional override value for the environment variable $RUN1D")
@click.option("--custom", type=str, help="Name of Custom Coadd schema")
@click.option("-p", "--plot", is_flag=True, help="Create spec plots")
@click.option("-e", "--epoch", is_flag=True, help="Run for epoch Coadds")
@click.option("--lsdr10", is_flag=True, help="Include Legacy Survey DR10 links on HTML")
@click.option("--allsky", is_flag=True, help="Reformat for Allsky Custom Coadd")
@click.pass_context
def run_reformat(ctx, **kwrds):
    """Build Spec Files"""
    args = AttrDict(ctx.params)
    spSpec_reformat(args.topdir, args.run2d, args.run1d, args.field, args.mjd,
                    plot=args.plot, epoch=args.epoch, lsdr10=args.lsdr10,
                    allsky=args.allsky, custom=args.custom)


@run.command(name='spcalib_qa', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.option('--run2d', envvar='RUN2D', help='Optional override value for the enviro variable $RUN2D')
@click.option('--bsr', 'boss_spectro_redux', envvar='BOSS_SPECTRO_REDUX', 
              help='Optional override value for the enviro variable $BOSS_SPECTRO_REDUX')
@click.option('--field', '-f',  type=str, help='Run for a single Field', default=None)
@click.option('--mjd',   '-m',  type=str, help='Run for a single MJD', default=None)
@click.option('--rerun', '-r',  is_flag = True, default=False, help='Rerun for all field-mjds in spAll')
@click.option('--nobkup', '-n', is_flag = True, default=False, help='Do not backup output and log file')
@click.option('--epoch',  '-e', is_flag = True, default=False, help='run for epoch coadds')
@click.option('--catchup','-c', is_flag = True, default=False, help='Run for missing field-mjds spAll')
@click.option('--outdir', help='Location to Save plots to (overrides the defaults)')
@click.option('--run2d_alt', 'run2d_alt', help='Alternative RUN2D for comparison')
@click.option('-bsra','--boss_spectro_redux_alt', 'boss_spectro_redux_alt', 
              help='Alternative BOSS_SPECTRO_REDUX for comparison')
@click.option('--plot_only', is_flag=True, default=False, help='Create Plots but dont update summary file')

@click.pass_context
def run_calibqa(ctx, **kwrds):
    """Compare photometric accuracy of standards"""
    args = AttrDict(ctx.params)
    spcalib_qa(**args)


@run.command(name='update_flags', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.option("--topdir", type=str, envvar='BOSS_SPECTRO_REDUX',
              show_default="env: BOSS_SPECTRO_REDUX",
              help="Optional override value for the environment variable $BOSS_SPECTRO_REDUX")
@click.option("--run2d", type=str, envvar="RUN2D", show_default="env: RUN2D",
              help="Optional override value for the environment variable $RUN2D")
@click.option("--custom", type=str, help="Name of Custom Coadd schema")
@click.option("--clobber", is_flag=True, help="Clobber spTargeting file")
@click.option("--nobackup", is_flag=True, help="Skip backup of existing summary files")
@click.pass_context
def run_updateflags(ctx, **kwrds):
    """Update SDSSV Targeting flats in the summary files"""
    args = AttrDict(ctx.params)
    if args.run2d is None:
        args.run2d = os.getenv('RUN2D')
    if args.topdir is None:
        args.topdir = os.getenv('BOSS_SPECTRO_REDUX')
    
    update_Targeting_flags(args.run2d, args.topdir, schema=args.custom, 
                           clobber=args.clobber, nobackup=args.nobackup)


@run.command(name='Plot_QA', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.option("-r", "--run2d", multiple=True, show_default="env: RUN2D", help="List of run2ds",
              default=lambda: (os.getenv("RUN2D"),) if os.getenv("RUN2D") else ())
@click.option("--test", multiple=True, callback=multi_str2bool,
              help="List of True/False test run2d (corresponding to run2d)")
@click.option("--test_path", default="/test/sean/", show_default=True,
              help="test Run2d path modification")
@click.option("--mjds_low", multiple=True, callback=multi_str2none,
              help="List of mjd lower limits (use 'None' for no limit)")
@click.option("--mjds_high", multiple=True, callback=multi_str2none,
              help="List of mjd upper limits (use 'None' for no limit)")
@click.option("--clobber_lists", is_flag=True, help="Clobber list of fieldIDs")
@click.option("--lco/--apo", "lco", default=False, help="Flag for LCO vs APO")
@click.option("--publish", is_flag=True, help="Create publication version of plot")
@click.option("--html", is_flag=True,
              help="Produces Plotly interactive HTML versions of the plots")
@click.option("--html_name", default=None,
              help="Name of HTML file (default = BOSS_QA-{obs}.html)")
@click.option("-f", "--fast_opsdb", is_flag=True,
              help="Skips OpsDB queries for SOS SN2 (and uses cached if available)")
@click.option("-e", "--epoch", is_flag=True, help="Produce plots for epoch coadds")
@click.option("-c", "--cron", is_flag=True, help="Produce cronlogs")
@click.option('--fid','--fieldids','fieldids', multiple=True, default=(), 
              help='Limit to these FieldIDs')
@click.option('--compare', is_flag = True, default=False, help="Direct Comparison of the run2ds in list")
@click.option('--start_mjd', type=int, default= None, help='Limit to MJDs on or after')
@click.option('--end_mjd', type=int, default= None, help='Limit to MJDs on or before')
@click.option('--output_dir', default= None, help='Overrides the output directory')

@click.pass_context
def run_plotqa(ctx, **kwrds):
    """Plot the SpectroPhotometry and SN2 QA plots"""
    args = AttrDict(ctx.params)
    if len(args.fieldid): args.fieldid = None
    mjds = {}

    args.run2d = list(args.run2d)
    args.test = list(args.test)
    args.mjds_low = list(args.mjds_low)
    args.mjds_high = list(args.mjds_high)

    n = len(args.run2d)
    args.test = args.test + [False] * (n - len(args.test))
    args.mjds_low = args.mjds_low + [None] * (n - len(args.mjds_low))
    args.mjds_high = args.mjds_high + [None] * (n - len(args.mjds_high))

    for i, run2d in enumerate(args.run2d):
        print(i,run2d)
        if args.mjds_high[i] == 'None':
            args.mjds_high[i] = None
        elif args.mjds_high[i] is None:
            args.mjds_high[i] = None
        else:
            args.mjds_high[i] = int(args.mjds_high[i])
        if args.mjds_low[i] == 'None':
            args.mjds_low[i] = None
        elif args.mjds_low[i] is None:
            args.mjds_low[i] = None
        else:
            args.mjds_low[i] = int(args.mjds_low[i])
        mjds[run2d] = [args.mjds_low[i], args.mjds_high[i]]
    obs='LCO' if args.lco is True else 'APO'
        

    plot_QA(args.run2d, args.test, mjds=mjds, obs=obs, testp=args.test_path,
            clobber_lists=args.clobber_lists, publish= args.publish,
            epoch=args.epoch, html_name=args.html_name, html=args.html,
            fast_opsdb = args.fast_opsdb, cron = args.cron)

if __name__ == "__main__":
    run()