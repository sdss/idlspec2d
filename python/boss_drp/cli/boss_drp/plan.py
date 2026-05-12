from sdss_access import Access
from boss_drp.prep.spplan import spplan2d, spplan1d
from boss_drp.prep.spplan_epoch import spplancombin
from boss_drp.prep.manage_coadd_Schema import manage_coadd_Schema
from boss_drp.prep.spplan_trace import spplanTrace
from boss_drp.prep.spplan_target import batch, CustomCoadd
from boss_drp.field.generations import generations
from boss_drp.Config import config, update_key, show_config, show_config_opt, fill_none_with_false
from boss_drp.utils.argparse_help import full_help_callback, OrderedGroup, _add_obs, AttrDict


import numpy as np
import click

@click.group('plan', cls = OrderedGroup,
             context_settings={"help_option_names": ['-h','--help'], 
                               "max_content_width": 150}) 
@click.option(
    "--fullhelp",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=full_help_callback,
    help="Show full help including all subcommands"
)
def plan():
    """BOSS DRP Plan Commands"""
    pass


def config_options(pipe_def="boss_drp"):
    def decorator(f):
        f = click.option( "--pipe_config", "--pc", "config", 
                         default=pipe_def, help="Queue Config name")(f)
        f = click.option("--pipe_config_file", "--pcf", "config_file",
                         default=None, help="Queue Config File Path")(f)
        f = show_config_opt(f)
        return f
    return decorator


def general_options(log_var = 'dailyplan_logfile', verbose_var='daily_plan_verbose', show_verbose=True, show1d=True, show_man=True):
    def decorator(f):
        f = click.option("--topdir", "BOSS_SPECTRO_REDUX",
                        help="Base run2d directory to BOSS_SPECTRO_REDUX environmental variable")(f)
        f = click.option("--run2d", "RUN2D",
                        help="Run2d to environmental variable")(f)
        if show1d:
            f = click.option("--run1d", "RUN1D",
                             help="Run1d to environmental variable")(f)
        # If you want repeated flags like --obs apo --obs lco
        f = click.option("--obs", "obs", multiple=True, default= (),
                        type=click.Choice(["apo", "lco"], case_sensitive=False),
                        help="Observatory {apo,lco}", )(f)

        # These are the Click-native equivalents of append_const to the same dest.
        f = click.option("--apo", is_flag=True, expose_value=False, help="Run apo",
                         callback=lambda ctx, param, value: _add_obs(ctx, param, value, "apo"))(f)

        f = click.option("--lco", is_flag=True, expose_value=False, help="Run lco",
                         callback=lambda ctx, param, value: _add_obs(ctx, param, value, "lco"))(f)

        f = click.option("--logfile", log_var, 
                        help="Optional logfile (Including path)")(f)
        if show_verbose:
            f = click.option("--verbose", verbose_var, 
                            help="Provide information about nonutlized frames")(f)
        f = click.option("-c", "--clobber", "clobber_plan",
                        is_flag=True, help="overwrites previous plan file")(f)
        if show_man:
            f = click.option("--override_manual/--no-override_manual", default=None,
                            help="Override/clobber manually edited plan")(f)
        return f
    return decorator

def sdss_access_options(f):
    f = click.option("--release", "RELEASE", default="sdsswork",
                     help="sdss_access data release")(f)
    f = click.option("--remote/--no-remote", "REMOTE", default=None,
                     help="allow for remote access to data using sdss-access")(f)
    return f
def Mos_targ_options(f):
    f = click.option("--v_targ", "V_TARG", default=None, #
                     help='SDSS-V MOS Targeting Product Version (for no Database access use)')(f)
    return f

def step_options(f):
    f = click.option("--skip2d/--no-skip2d", default=None, help="Skip spplan2d")(f)
    f = click.option("--skip1d/--no-skip1d", default=None, help="Skip spplan1d")(f)
    return f

def filter_options(mjd=False, field=False, generation =False, devel=False, trace=False):
    def decorator(f):
        if mjd:
            f = click.option("--mjd", multiple=True, help="Use data from these MJDs.")(f)
            f = click.option("--mjdstart", help="Starting MJD")(f)
            f = click.option("--mjdend", help="Ending MJD")(f)

        if field:
            f = click.option("--field", multiple=True, help="Use data from these fields.")(f)
            f = click.option("--fieldstart", help="Starting Field")(f)
            f = click.option("--fieldend", help="Ending Field")(f)

        if generation:
            f = click.option("--legacy/--no-legacy", default=None, help="Include legacy (BOSS/eBOSS) plates")(f)
            f = click.option("--plates/--no-plates", default=None, help="Include SDSS-V plates")(f)
            f = click.option("--fps/--no-fps", default=None, help="Include FPS Fields")(f)
            f = click.option("--sdssv/--no-sdssv", default=None, help="Include both SDSS-V Fields & Plates")(f)
        if devel:
            f = click.option("--commissioning/--no-commissioning", default=None, help="Include SDSS-V FPS Commission Fields")(f)
            f = click.option("--dither/--no-dither", "dither", default=None, help="Include Dither fields")(f)

        if trace:
            f = click.option('--mjd_plans/--no-mjd_plans', 'trace_all_mjds', default=None, help='Only build plans for MJDs with spPlan2d')(f)
            f = click.option('--include_hartmann', default=False, is_flag=True, 
                             help=('Include Hartmann exposures in spPlanTrace as Arc frames (only relevent for testing and not '+
                                   'recommended for production since Hartmanns are not real Arc exposures)'))(f)
            f = click.option('--exclude_arc', default=False, is_flag=True, 
                        help=('Exclude Arc exposures in spPlanTrace (only relevent for testing and not '+
                            'recommended for production since Arcs are needed)'))(f)
        return f
    return decorator

def run2d_options(f):
    f = click.option("--matched_flats/--no-matched_flats", default=None, help="Require Flat from a field/plate")(f)
    f = click.option("--matched_arcs/--no-matched_arcs", default=None, help="Allow Arc from another field/plate")(f)
    f = click.option("--minexp", "minscience", default=1, type=int,
                     help="Min Science Exposures in Plan (default=1)")(f)
    f = click.option("--multiple_flat/--no-multiple_flat", default=None, help="Find all possible flat calibration frames")(f)
    f = click.option("--multiple_arc/--no-multiple_arc", default=None, help="Find all possible arc calibration frames")(f)
    f = click.option("--manual_noarc/--no-manual_noarc", "flag_nomatch_manual", default=None,
                     help="if nomatched_arcs is False, builds spplan with unmatched arcs and mark as manual")(f)
    return f


def run1d_options(f):
    f = click.option("--plate_epoch/--no-plate_epoch", default=None,
                     help="Use a variable max epoch length for plate coadd")(f)
    f = click.option("--quick/--no-quick", "quick1d", default=None,
                     help="Use the list of new spPlan2d as a filter for fields")(f)
    return f

def epoch_options(f):
    f = click.option("--minexp", "minscience", default=1, type=int,
                     help="Min Science Exposures in Plan (default=1)")(f)
    f = click.option('--abandoned/--no-abandoned', default=None, help="Create plans for abandoned epochs")(f)
    f = click.option('--started/--no-started', default=None, help="Create plans for started epochs (including unfinished)")(f)
    f = click.option('--min_epoch_len', help="minimum length of epoch required to produce plan", type=int, default = None)(f)
    return f

def target_options(f):
    f = click.option("--batch/--manual", "batch", default=None,
                      help="Batch run all active Coadd Schema (default) or manually run a single schema")(f)

    f = click.option('--name', 'custom_name', help = 'Name of Custom Coadd')(f)
    f = click.option('--coaddfile', 'schema_file', help = 'File of store Coadding Schema')(f)

    f = click.option('--DR/--no-DR', default=None, help = 'DR/IPL Batch Coadding')(f)
    f = click.option('--cartons', multiple=True, help = "list of cartons")(f)
    f = click.option('--catalogids', multiple=True, help = "list of sdss_ids (or catalogids)")(f)
    f = click.option('--program', multiple=True, help = 'list of programs')(f)
    f = click.option('--coadd_mjdstart', help = 'First Coadd MJD to include')(f)
    f = click.option('--rerun1d/--no-rerun1d', default=None,
                     help = 'Provides flag for coadd to be rerun though 1D analysis')(f)
    f = click.option('--use_catid/--no-use_catid', '-u', default=None,
                      help='Uses CatalogID rather then sdss_id')(f)
    f = click.option('--use_firstcarton/--no-use_firstcarton', default=False, 
                     help='Use Firstcarton only for carton match (dont look at db)')(f)
    f = click.option('--useDB/--no-useDB', default=False, 
                     help='Use sdss targetdb instead of the Semaphore targeting flag (if not use_firstcarton)')(f)
    return f


def check_release(RELEASE, REMOTE):
    if RELEASE != 'sdsswork':
        if RELEASE not in Access().get_available_releases():
            raise click.ClickException(f"{RELEASE} is not a valid release")
    else:
        if REMOTE is True:
            try:
                Access().remote()
            except:
                raise click.ClickException('ERROR: No netrc file found. see https://sdss-access.readthedocs.io/en/latest/auth.html#auth')
    return

def load_config(exclude_args=None, warn=True, **args):
    config.load(config_name=args.get('config',None), config_file=args.get('config_file',None))

    if exclude_args is None:
        exclude_args = []
    exclude_args.extend(['config','config_file','show_config'])
    for arg, value in args.items():
        if value is not None and arg not in exclude_args:
            if not update_key(config.pipe, arg, value):
                if warn:
                    print(f"[DEBUG] Key '{arg}' not found anywhere in config.")
    
    if config.pipe['general.V_TARG'] is None:
        update_key(config.pipe, 'V_TARG', '*')
        
    return



@plan.command(name='daily', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@config_options(pipe_def = 'boss_drp')
@step_options
@general_options(show1d=False)
@sdss_access_options
@Mos_targ_options
@filter_options(mjd=True, field=True, generation=True, devel=True)
@run2d_options
@run1d_options
@click.pass_context
def daily(ctx, **kwrds):
    """Produce the spPlan2d and spPlancomb files for the pipeline run"""
    args = AttrDict(ctx.params)
    
    load_config(**args)

    obs = config.pipe['fmjdselect.obs'] or ['apo','lco']
    obs = np.atleast_1d(obs).tolist()
    generations.set_config(obs[0])
    for i,_obs in enumerate(obs):
        if i == 1:
            generations.set_config(_obs)
            load_config(warn=False,**args)

        update_key(config.pipe, 'obs', _obs)
        fill_none_with_false(config.pipe)
        if args.show_config:
            show_config(queue=False)
            continue
        plans2d = None
        
        if not config.pipe['plan.daily.skip2d']:
            plans2d = spplan2d()
            if config.pipe['plan.daily.quick1d']:
                if plans2d is None:
                    update_key(config.pipe, 'skip1d', True)
        if not config.pipe['plan.daily.skip1d']:
            spplan1d(plans = plans2d)

@plan.command(name='trace', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@config_options(pipe_def = 'boss_drp')
@general_options(verbose_var='traceplan_verbose', log_var='traceplan_logfile', show1d=False)
@sdss_access_options
@filter_options(mjd=True, trace=True)
@click.pass_context
def trace(ctx, **kwrds):
    """Produces spPlanTrace for the Use of Master Arc and Flat Frames to build Traces"""
    args = AttrDict(ctx.params)
    args['clobber_spTrace'] = args.clobber_plan
    check_release(args['RELEASE'], args['REMOTE'])
    exclude=['include_hartmann', 'exclude_arc']
    load_config(**args, exclude_args=exclude)
    
    obs = config.pipe['fmjdselect.obs']
    obs = np.atleast_1d(obs).tolist()
    generations.set_config(obs[0])

    for i, _obs in enumerate(obs):
        if i == 1:
            generations.set_config(_obs)
            load_config(warn=False,**args)

        update_key(config.pipe, 'obs', _obs)
        fill_none_with_false(config.pipe)
        if args.show_config:
            show_config(queue=False)
            continue
        spplanTrace(include_hartmann=args.include_hartmann, exclude_arc = args.exclude_arc) 

@plan.command(name='epoch', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@config_options(pipe_def = 'boss_drp')
@general_options(show_verbose=False, log_var='epochplan_logfile')
@sdss_access_options
@Mos_targ_options 
@filter_options(mjd=True, field=True, generation=True)
@epoch_options
@click.pass_context
def epoch(ctx, **kwrds):
    """Builds the spPlancombepoch files for the Epoch Coadd Pipeline Runs"""
    args = AttrDict(ctx.params)
    
    check_release(args['RELEASE'], args['REMOTE'])
    load_config(**args)
    
    obs = config.pipe['fmjdselect.obs']
    obs = np.atleast_1d(obs).tolist()
    generations.set_config(obs[0])

    for i, _obs in enumerate(obs):
        if i == 1:
            generations.set_config(_obs)
            load_config(warn=False,**args)
    
        update_key(config.pipe, 'obs', _obs)
        fill_none_with_false(config.pipe)
        if args.show_config:
            show_config(queue=False)
            continue
        spplancombin()

@plan.command(name='CoaddSchema', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.option("--coaddfile", "-f", default=None,
              help="File to store Coadding Schema (Default: {topdir}/{run2d}/fields/SDSSV_BHM_COADDS.par)")
@click.option("--topdir", envvar='BOSS_SPECTRO_REDUX',help="Override value for the environment variable $BOSS_SPECTRO_REDUX.")
@click.option("--run2d", envvar='RUN2D', help="Override value for the environment variable $RUN2D")
@click.option("--name", default=None, help="Name of Custom Coadd")
@click.option("--DR", is_flag=True, help="DR/IPL Coadding")
@click.option("--rerun1d", "-r", is_flag=True, help="Provides flag for coadd to be rerun though 1D analysis")
@click.option("--active", "-a", is_flag=True, help="Activate (or deactivate) a Coadding Schema")
@click.option("--carton", "-c", multiple=True, help="list of cartons")
@click.option("--SDSSIDS", "-i", multiple=True, help="list of SDSS_IDS (or CatalogIDs if use_catid is set)")
@click.option("--program", "-p", multiple=True, help="list of programs")
@click.option("--legacy", "-l", multiple=True, help="list of Legacy Tags to include")
@click.option("--use_catid", "-u", is_flag=True, help="Use CatalogIDs rather then SDSS_IDs")
@click.option("--use_firstcarton", is_flag=True, help="Use Firstcarton only for carton match (dont look at db)")
@click.option("--cadence", "-t", type=float, default=0.0, help="Number of days between coadd epochs")
@click.option("--show", "-s", is_flag=True, help="Show Configurations")
@click.option("--mjd", multiple=True, help="Use data from these MJDs.")
def run_CoaddSchema(coaddfile, topdir ,run2d, name, DR, rerun1d, active, carton, SDSSIDS,
                    program, legacy, use_catid, use_firstcarton, cadence, show, mjd):
    """Manage SDSSID/Catalogid Custom Coadds Schema"""
    manage_coadd_Schema(name, topdir=topdir, run2d=run2d, DR=DR, CARTON=carton, CATID=SDSSIDS,
                        PROGRAM=program, RERUN1D=rerun1d, CADENCE=cadence, MJD=mjd, ACTIVE=active, 
                        legacy = legacy, coaddfile=coaddfile, show=show, use_catid=use_catid,
                        use_firstcarton=use_firstcarton)


@plan.command(name='target', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@config_options(pipe_def = 'boss_drp')
@general_options(log_var='customplan_logfile', show_verbose=False, show_man=False) 
@filter_options(mjd=True)
@click.pass_context
def target(ctx, **kwrds):
    """Build SDSSID/CatalogID Custom Combine Plan"""
    args = AttrDict(ctx.params)
    check_release(args['RELEASE'], args['REMOTE'])
    load_config(**args)
    generations.set_config(config.pipe['fmjdselect.obs'])

    if ((not config.pipe['plan.custom.DR']) and (config.pipe['plan.custom.cartons'] is None) and
        (config.pipe['fmjdselect.mjd'] is None) and (config.pipe['fmjdselect.mjdstart'] is None) and
        (config.pipe['fmjdselect.mjdend'] is None) and (not config.pipe['plan.custom.rerun1d']) and
        (config.pipe['plan.custom.program'] is None)):

        update_key(config.pipe,'batch', True)
    fill_none_with_false(config.pipe)
    if args.show_config:
        show_config(queue=False)
        return
    if config.pipe['plan.custom.batch']: 
        batch()
    else:
        CustomCoadd(config.pipe['customSettings.custom_name'],config.pipe['general.BOSS_SPECTRO_REDUX'],
                    config.pipe['general.RUN2D'],config.pipe['general.RUN1D'], cartons = config.pipe['plan.custom.cartons'],
                    catalogids = config.pipe['plan.custom.catalogids'], obs = config.pipe['fmjdselect.obs'], 
                    clobber = config.pipe['Clobber.clobber_plan'], logfile = config.pipe['plan.custom.customplan_logfile'],
                    mjd = config.pipe['fmjdselect.mjd'], mjdstart = config.pipe['fmjdselect.mjdstart'],
                    mjdend = config.pipe['fmjdselect.mjdend'], program = config.pipe['plan.custom.program'],
                    rerun1d = config.pipe['plan.custom.rerun1d'], use_catid = config.pipe['plan.custom.use_catid'],
                    use_firstcarton=config.pipe['plan.custom.use_firstcarton'],
                    coadd_mjdstart=config.pipe['plan.custom.coadd_mjdstart'], useDB=config.pipe['plan.custom.useDB'])
