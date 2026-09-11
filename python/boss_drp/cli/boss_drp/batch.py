import click

### These are now imported lazily in the functions that need them
#from boss_drp.run.slurm_readfibermap import slurm_readfibermap
#from boss_drp.run.slurm_runfix import slurm_runfix
#from boss_drp.run.slurm_sos import slurm_SOS
#from boss_drp.run.slurm_Summary import slurm_Summary
#from boss_drp.run.slurm_spTrace import run_spTrace as slurm_run_spTrace
#from boss_drp.run.uubatchpbs import uubatchpbs
########

from boss_drp.Config import config, show_config_opt, fill_none_with_false, show_config, update_key
from boss_drp.utils.argparse_help import full_help_callback, AttrDict, _add_obs
from boss_drp.cli.boss_drp.cli2config import cli2config
from boss_drp.utils import jdate
import os

@click.group(name='batch',context_settings={"help_option_names": ['-h','--help']})
@click.option(
    "--fullhelp",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=full_help_callback,
    help="Show full help including all subcommands"
)
def batch():
    """BOSS DRP Batch Cluster Tools"""
    pass


def config_opts(default_qc = 'readfibermap', sos=False):
    def decorator(f):
        if not sos:
            f = click.option('--pipe_config', '--pc', 'config', default='boss_drp', help='Queue Config name')(f)
            f = click.option('--pipe_config_file', '--pcf', 'config_file', default=None, help='Queue Config File Path')(f)
        f = click.option('--queue_config', '--qc', default=default_qc, help='Queue Config name')(f)
        f = click.option('--queue_config_file', '--qcf', default=None, help='Queue Config File Path')(f)
        f = show_config_opt(f)
        return f
    return decorator

def obs_opts(f):
    f = click.option('--obs', 'obs', multiple=True, type=click.Choice(['apo', 'lco'], case_sensitive=False),
                     help='Observatory {apo,lco}. Can be repeated.')(f)
    f = click.option("--apo", is_flag=True, expose_value=False, help="Run apo",
                        callback=lambda ctx, param, value: _add_obs(ctx, param, value, "apo"))(f)
    f = click.option("--lco", is_flag=True, expose_value=False, help="Run lco",
                        callback=lambda ctx, param, value: _add_obs(ctx, param, value, "lco"))(f)
    return f

def mjd_opts(daily=False):
    def decorator(f):
        f = click.option('--mjd', multiple=True, type=int, help='MJD dates to reduce; default="*"')(f)
        f = click.option('--mjdstart', type=int, default=None, help='Starting MJD')(f)
        if daily:
            f = click.option('--daily/--no-daily', 'trace_all_mjds', is_flag=True, default=None,
                            help='Run in daily mode (only use MJDs specified)')(f)
        f = click.option('--mjdend', type=int, default=None, help='Ending MJD')(f)
        return f
    return decorator

def queue_opts(maxjobs=False):
    def decorator(f):
        f = click.option('--mem_per_cpu', type=str, default=None, help='Memory allocated per CPU')(f)
        f = click.option('--walltime', type=str, default=None, help='Wall time in hours')(f)
        f = click.option('--nodes', type=int, default=None, help='Number of Nodes')(f)
        if not maxjobs:
            f = click.option('--ppn', type=int, help='Number of processors per node')(f)
        else:
            f = click.option('--maxjobs', type=int, help='Max Number of Parallel jobs per node')(f)
        f = click.option('--nbundle', type=int, default=None, help='Number of jobs to bundle')(f)
        return f
    return decorator

@batch.command(name='readfibermap',context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
# Basic options
@click.option('--topdir', default=None, help='Boss Spectro Redux base directory')
@click.option('--run2d', default=None, help='Run2d')
@click.option('--clobber/--no-clobber', 'clobber_fibermap', is_flag=True, default=None, help='Clobber spfibermaps')
@obs_opts
@click.option('--v_targ', '--V_TARG', 'V_TARG',
              help='SDSS-V MOS Targeting Product Version (for no Database access use)')
@mjd_opts()
@config_opts()
@queue_opts()
@click.pass_context
def run_readfibermap(ctx, **kwrds):
    """Create a batch readfibermap job. Without access to the SDSS Slurm package, it prints the commands for manual execution"""
    from boss_drp.run.slurm_readfibermap import slurm_readfibermap
    args = AttrDict(ctx.params)

    if len(args.obs)== 0:
        args.obs = ['apo','lco']

    config_par = {'mem_per_cpu':args.mem_per_cpu,
                  'wall':args.walltime,
                  'ppn': (args.ppn or 20),
                  'nbundle':args.nbundle}

    cli2config(args, config_par = config_par)
    fill_none_with_false(config.pipe)
    fill_none_with_false(config.queue)
    if args.show_config:
        show_config()
        return
    slurm_readfibermap()



@batch.command(name='runfix', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.option('--full/--no-full', is_flag=True, default=None, help='Rerun full pipeline regardless of crashed step')
@click.option('--running/--no-running', is_flag=True, default=None, help='Select Field-MJDs with running status')
@click.option('--topdir', 'BOSS_SPECTRO_REDUX', type=str, default=None,
              help='Optional override value for the config')
@click.option('--run2d', 'RUN2D', type=str, default=None,
              help='Optional override value for the config')
@click.option('--run1d', 'RUN1D', type=str, default=None, 
              help='Optional override value for the config')
@click.option('--epoch/--daily', is_flag=True, default=None, help='Run for epoch Coadds')
@obs_opts
@mjd_opts()
@config_opts(default_qc='tagged_daily')
@queue_opts()
@click.option('--no-write/--write','no_write', is_flag=True, default=None, help='skip writing and submitting job')
@click.option('--nosubmit/--submit', 'no_submit',is_flag=True, default=None, help='Build, but not submit redux files')
@click.pass_context
def run_runfix(ctx, **kwrds):
    """Check for failed runs and setup the runs to clean and rerun the crashed field-mjds"""
    from boss_drp.run.slurm_runfix import slurm_runfix
    args = AttrDict(ctx.params)
    args.custom = None
    
    if not args.obs:
        args.obs = ['apo','lco']
    config_par = {'no_write': args.no_write,
                  'wall':args.walltime,
                  'nodes':args.nodes,
                  'no_submit': args.no_submit,
                  'mem_per_cpu':args.mem_per_cpu,
                  'nbundle':args.nbundle}
    cli2config(args, config_par = config_par, exclude=['full','running','no_submit'])
    fill_none_with_false(config.pipe)
    fill_none_with_false(config.queue)

    if args.show_config:
        show_config()
        return
    slurm_runfix(full=args.full, fix_running=args.running)

@batch.command(name='sos', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@obs_opts
@mjd_opts()
@click.option('--no-reject','no_reject', is_flag=True, default=False,
              help="Overrides the Calibration rejection (use with caution)")
@click.option('--clobber_fibermap', '-f', is_flag=True, default=False,
              help="Clobbers the existing spfibermap files")
@click.option('-n', '--no-arc2trace','no_arc2trace', is_flag=True, default=False,
             help="Skip Utilizing arc2trace refinements")
@click.option('-o', '--forcea2t', is_flag=True, default=False,
              help="Force arc2trace for all fields (even if flat exists for field)")
# Boolean pair (store_false / store_true)
@click.option('--sdssv-sn2/--no-sdssv-sn2', 'sdssv_sn2', default=True,
              help="Report a second set of SN2 values with updated fit parameters")
@click.option('--sn2-15/--no-sn2-15', 'sn2_15', default=True,
              help="Report a set of SN2 values with a fiducial mag of 15")
@click.option('--bright', is_flag=True, default=True, help='Display BOSS_only Bright Time Operation SN2_15')
@config_opts(default_qc='sos', sos=True)
@queue_opts()
@click.option('--no-submit', is_flag=True, default=False, help='Build, but not submit redux files')
@click.pass_context
def run_sos(ctx, **kwrds):
    """Create SOS queue job. Without access to the SDSS Slurm package, it prints the commands for manual execution"""
    from boss_drp.run.slurm_sos import slurm_SOS
    args = AttrDict(ctx.params)
    obs = args.obs
    args.config = None
    args.config_file = None
    config_par = {'wall':args.walltime,
                  'ppn': args.ppn,
                  'nodes': args.nodes,
                  'mem_per_cpu': args.mem_per_cpu,
                  'no_submit': args.no_submit,
                  'nbundle':args.nbundle}
    args.queue_config = None
    args.queue_config_file = None
    cli2config(args, config_par = config_par)
    fill_none_with_false(config.pipe)
    fill_none_with_false(config.queue)

    args.obs = obs
    if args.show_config:
        show_config(pipe=False)
        return
    #TODO - use config????
    slurm_SOS(**args)


@batch.command(name='spTrace', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150),
               short_help='Create spTrace Queue jobs') 
@click.option('--topdir', 'BOSS_SPECTRO_REDUX', type=str, default=None,
              help='Optional override value for the config')
@click.option('--run2d', 'RUN2D', type=str, default=None,
              help='Optional override value for the config')
@obs_opts
@click.option('--clobber/--no-clobber', 'clobber_spTrace', is_flag = True, default=None, 
              help='Clobber the existing Plan files')
@click.option('--debug/--no-debug', 'debug', is_flag = True, default=None, 
              help='Save sdssproc outputs')
@click.option('--saveraw/--no-saveraw', 'saveraw', is_flag = True, default=None, 
              help='Clobber the existing Plan files')
@click.option('--no-skip_plan/--skip_plan', 'run_spTrace_plan', is_flag = True, default=None, 
              help='Skip creating plans and use currently existing plans')
@click.option('--hartmann', '--hart','hartmann', is_flag=True, default=False, help='Trace and Extract the Hartmann Frames')
@mjd_opts(daily=True)
@config_opts(default_qc='spTrace')
@queue_opts(maxjobs=True)
@click.option('--nosubmit/--submit', 'no_submit', is_flag=True, default=None, help='Build, but not submit redux files')
@click.pass_context
def run_spTrace(ctx, obs, **kwrds):
    """
    Create spTrace Queue jobs. Without access to the SDSS Slurm package, it prints the commands for manual execution.
    """
    from boss_drp.run.slurm_spTrace import run_spTrace as slurm_run_spTrace
    args = AttrDict(ctx.params)
 
    
    if not args.obs:
        args.obs = ['apo','lco']


    if args.debug:
        args.saveraw= True

    runobs = args.obs
    for obs in runobs:
        args.obs = obs

        config_par = {'wall':args.walltime,
                    'nodes':args.nodes,
                    'no_submit':args.no_submit,
                    'max.jobs':args.maxjobs,
                    'mem_per_cpu':args.mem_per_cpu,
                    'nbundle':args.nbundle}

        cli2config(args, config_par = config_par, exclude=['hartmann'])

        if config.pipe['fmjdselect.mjd'] is None:
            if config.pipe['fmjdselect.mjdrange'] is None:
                update_key(config.pipe,'mjdrange',[jdate.astype(int),jdate.astype(int)+1])

        fill_none_with_false(config.pipe)
        fill_none_with_false(config.queue)

        if args.show_config:
            show_config()
            continue

        slurm_run_spTrace(hartmann=args.hartmann)


@batch.command(name='Summary', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.option('--module', '-m', default = None, help = 'module file to use (ex bhm/master[default] or bhm/v6_0_9)')
@click.option('--topdir', 'BOSS_SPECTRO_REDUX', type=str, default=None,
              help='Optional override value for the config')
@click.option('--run2d', 'RUN2D', type=str, default=None,
              help='Optional override value for the config')
@click.option('--run1d', 'RUN1D', type=str,default=None,
              help='Optional override value for the config')
@click.option('--epoch/--no-epoch', is_flag=True, default=None, help='Run for epoch Coadds')
@click.option('--allsky/--no-allsky', is_flag=True, default=None, help='Run for custom allsky Coadds')
@click.option('--custom', "custom_name", default=None, help='Run for epoch Coadds')
@click.option("--to_fits", is_flag=True, help="Dump Parquet to fits format")
@click.option("--keep_active", is_flag=True, help='Run "touch" on all intermediate files to keep them active')
@click.option("--force", "--force_rebuild", "force_rebuild", is_flag=True,
              help="Rebuild Summary even if nothing changed")
@click.option("--clobber_mjd", "clobber_mjd", is_flag=True, help="Clobber all spAll-MJD files")

@click.option('--daily/--no-daily', 'after_daily', is_flag=True, default=None, 
              help='only run if daily run has been run today')
@click.option('--monitor/--no-monitor', 'pipe_monitor', is_flag=True, default=None, 
              help='Monitor job and send email at completion with the logs')
@click.option('--fieldlist/--no-fieldlist', 'run_fieldlist', is_flag=True, default=None, 
              help='Running Fieldlist')
@click.option('--backup', 'backup', type=int, default=None,
              help='Number of backups to keep, or None (or 0) to not create backup')
@click.option('--n_iter', 'n_iter', type=int, default=None,
              help='number of iterations of field merge to run')
@click.option('--ndays', 'ndays', type=int, default=None,
              help='Limit spAll update to last ndays')
@click.option('--skip_specprimary', 'skip_specprimary', flag_value='skip',
              default=None, help='Skip calculation of Specprimary')
@click.option('--update_specprimary', 'skip_specprimary', flag_value='update',
              help='Only update new Specprimary')
@click.option("--update_target_flags/--no-update_target_flags", 
              "--tf/--no-tf", "update_target_flags", default=False,
              help="Use the spTargeting file to update the summary file to the latest Targeting Flags")
@click.option('--utah_daily/--no_utah', "database", is_flag=True, default=None, 
              help='Load tagged daily run into Pipelines.boss_drp database table')
@click.option('--verbose/--no-verbose', "verbose",  is_flag=True, default=None, 
              help='Run Fieldmerge with verbose')
@click.option('--email_start/--no-email_start', "email_start",  is_flag=True, default=None, 
              help='Send email at start of run')
@click.option('--defaults', "defaults",  is_flag=True, default=None, 
              help='Sets --merge_only  --backup 3  --monitor --update_specprimary --ndays 10 --qc summary_full')
@config_opts(default_qc='summary')

@click.option('--walltime', type=str, default=None, help='Wall time in hours')
@click.option('--mem', default=None, help = 'memory in bytes')
@click.option('--ppn', help='Number of processors per node', type=int)
@click.option('--nosubmit/--submit', "no_submit", is_flag=True, default=None, help='Create queue job but do not submit it')
@click.pass_context
def run_summary(ctx, **kwrds):
    """Create daily field merge queue job"""
    from boss_drp.run.slurm_Summary import slurm_Summary
    args = AttrDict(ctx.params)

    if args.defaults:
        if args.queue_config is None:
            args.queue_config = 'summary_full'
        args.merge_only = True
        if args.backup is None:
            args.backup = 3
        args.pipe_monitor = True
        args.update_specprimary = True
        args.update_target_flags = True                 
        if args.ndays is None:
            args.ndays = 10
    elif args.queue_config is None: 
        args.queue_config = 'summary'

    #TODO: Decide how to handle defaults with the config file. 
    # If the user specifies a queue_config, should that override the defaults set by --defaults? 
    # For now, if --defaults is set, it will override any conflicting values in the queue_config, 
    # but if --defaults is not set, then the queue_config will determine all values including defaults. 
    # This allows for flexibility in using pre-set queue configs while still allowing for quick overrides with --defaults.
    
    #TODO: make sure all updated fieldmerge flags are included here

    if args.backup == 0:
        args.backup = None


    config_par = {'mem':args.mem,
                  'wall':args.walltime,
                  'ppn': args.ppn,
                  'no_submit': args.no_submit}

    cli2config(args, config_par = config_par, set_gen=False, exclude=['walltime'])
    fill_none_with_false(config.pipe)
    fill_none_with_false(config.queue)

    if config.pipe['general.module'] is None:
        module = os.getenv('MODULE', default=None)
        if module is None:
            module = os.getenv('RUN2D', default=None)
            if module is None:
                module = 'bhm/master'
        update_key(config.pipe, 'module', module)

    if args.show_config:
        show_config()
        return

    slurm_Summary()

def pipeline_options(f):
    # Short cuts
    f = click.option('--sdssv', is_flag=True, default=False,
        help='--mwm --no-reject --merge3d')(f)
    f = obs_opts(f)
    f = click.option('--bay15', 'map3d', flag_value='bayestar15', help='Set map3d to bayestar15 model')(f)
    # f = click.option('--eden23', 'map3d', flag_value='edenhofer2023', help='Set map3d to edenhofer2023 model')(f)
    f = click.option('--merge3d', 'map3d', flag_value='merge3d', help='Set map3d to best 3d model')(f)
    f = click.option('--topdir', 'BOSS_SPECTRO_REDUX', type=str, default=None,
                     help='Optional override value for the config')(f)
    f = click.option('--run1d', 'RUN1D', type=str, default=None,
                     help='Optional override value for the config')(f)
    f = click.option('--run2d', 'RUN2D', type=str, default=None,
                     help='Optional override value for the config')(f)
    f = click.option('--idlutils_1d', type=str, default=None,
                     help='idlutils override version of spec1d')(f)
    f = click.option('--no-reject/--reject','no_reject', is_flag=True, default=None,
                     help='Deactivate Rejection in Coadd')(f)
    f = click.option('--MWM-fluxer/--no-MWM-fluxer', '--mwm/--no-mwm', 'MWM_fluxer', is_flag=True, default=None, help='')(f)
#    f = click.option('--map3d', type=click.Choice(['bayestar15', 'bay15', 'merge3d', 'eden23'], case_sensitive=False),
    f = click.option('--map3d', type=click.Choice(['bayestar15', 'bay15', 'merge3d'], case_sensitive=False),
                    default=None, help='Name of 3d dustmap to use with MWM_fluxer (default=None)')(f)
    f = click.option('--noxcsao/--xcsao', 'run_XCSAO', is_flag=True, flag_value=False, default=None, help='Skip pyXCSAO')(f)
    f = click.option('--nodist/--dist', 'nodist', is_flag=True, default=None,
                     help='Unset --nodist and reactivate the flux distortion corrections')(f)
    f = click.option('--skip_specprimary','skip_specprimary', flag_value='skip', default=None,
                     help='Skip Calculation of Specprimary')(f)
    f = click.option('--update_specprimary','skip_specprimary', flag_value='update',
                     help='Only update new Specprimary')(f)
    f = click.option('--onestep_coadd/--no-onestep_coadd', is_flag=True, default=None,
                     help='Use legacy one step version of coadd')(f)
    f = click.option('--fibermap_clobber/--no-fibermap_clobber', 'clobber_fibermap', is_flag=True, default=None,
                     help='Clobber spfibermap fits file')(f)
    f = click.option('--saveraw/--no-saveraw', is_flag=True, default=None, help='Save sdssproc outputs')(f)
    f = click.option('--debug/--no-debug', is_flag=True, default=None, help='Save extraction debug files')(f)
    f = click.option('--no-db/--db', 'no_db', is_flag=True, default=None, help='skip Database operations')(f)
    f = click.option('--fast_no_db', required=False,
                     help='When using --no-db, streamlines process and only gets parallax from MOS target files')(f)
    f = click.option('--release', default=None, required=False, help='sdss_access data release ...')(f)
    f = click.option('--v_targ', '--V_TARG','V_TARG', default=None,
                     help='SDSS-V MOS Targeting Product Version  (for no Database access use)')(f)
    f = click.option('--a2t/--no-a2t', 'force_arc2trace', is_flag=True, default=None,
                     help='Force Use of Arc2Trace')(f)
    f = click.option('--clobber/--no-clobber', 'clobber_pipe', is_flag=True, default=None,
                     help='Clobber redux')(f)
    # Pipeline step options
    f = click.option('--no-healpix/--healpix','--nohp/--hp', 'run_healpix', is_flag = True, default=None, help='Turn off copy to healpix')(f)
    # f = click.option('--no-merge-spall/--merge_spall', 'run_Summarymerge', default=True, help='Skip building full SpAll File')(f)

    f = click.option('--skip2d/--no-skip2d', is_flag=True, default=None, help='Skip spreduce2d')(f)

    f = click.option('--only1d/--no-only1d', is_flag=True, default=None, help='Run spec1d step only (eg. spreduce1d_empca, XCSAO)')(f)

    # Select fields
    f = click.option('--field', '-f', multiple=True, type=str, help='Plate/Field numbers to reduce default="*"')(f)
    f = click.option('--fieldstart', default=None, type=str, help='Starting Field/Plate number')(f)
    f = click.option('--fieldend', default=None, type=str, help='End Field/Plate number')(f)

    # Select MJDs
    f = mjd_opts()(f)

    f = config_opts(default_qc = 'tagged_daily')(f)
    # Configurations

    # Queue options
    f = click.option('--no-write/--write', 'no_write', is_flag=True, default=None, help='skip writing and submitting job')(f)
    f = click.option('--mem-per-cpu','mem_per_cpu', type=str, help='Memory allocated per CPU')(f)
    f = click.option('--walltime', type=str, help='Wall time in hours')(f)
    f = click.option('--nodes', default=None, type=int, help='Number of Nodes')(f)
    f = click.option('--ppn', type=int, help='Number of processors per node')(f)
    f = click.option('--nosubmit/--submit', "no_submit", is_flag=True, default=None, help='Build, but not submit redux files')(f)
    f = click.option('--nbundle', type=int, default=None, help='Number of jobs to bundle')(f)

    # Custom coadd options
    f = click.option('--epoch/--no-epoch', is_flag=True, default=None, help='Epoch Coadds')(f)
    f = click.option('--custom', 'custom_name', type=str, help='Name of custom Coadd Schema')(f)
    f = click.option('--allsky/--no-allsky', is_flag=True, default=None, help='All Sky Coadds')(f)
    f = click.option('--single-mjd/--no-single-mjd', 'custom_single_mjd', is_flag=True, default=None,
                     help='Run Each Custom MJD coadd+1dpost as seperate job')(f)
    f = click.option('--coadd-only/--no-coadd-only','coadd_only', is_flag=True, default=None, help='Run spspec_target_merge only')(f)
    f = click.option('--1dpost/--all', 'post1d', is_flag=True, default=None, help='Run 1d analysis and post processing only')(f)

    # Email outputs
    f = click.option('--email/--no-email', 'sendemail', is_flag=True, default=None, help='Email log using $DAILY_DIR/etc/emails')(f)

    f = click.option('--allemail/--primary-email', is_flag=True, default=None,
                     help='Email intermediate log using all emails in $DAILY_DIR/etc/emails (defaults to first email only)')(f)

    return f

@batch.command(name='pipe', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150),
               short_help='Build idlspec2d redux and submit to the cluster queue.') 
@pipeline_options
@click.pass_context
def run_pipe(ctx, **kwrds):
    """
    Build idlspec2d redux and submit to the cluster queue. 
    Without access to the SDSS Slurm package, it prints the commands for manual execution
    """
    from boss_drp.run.uubatchpbs import uubatchpbs
    args = AttrDict(ctx.params)

    args.run_Summarymerge = False # This option is deprecated due to run time, but left here incase we ever want to add it back.
                                  # Currently the uubatchpbs code does not have this syntax in the template/redux Script.

    
    config_par = {'no_write': args.no_write,
                  'mem_per_cpu':args.mem_per_cpu,
                  'wall':args.walltime,
                  'nodes':args.nodes,
                  'ppn':args.ppn,
                  'no_submit': args.no_submit,
                  'nbundle':args.nbundle}

    cli2config(args, config_par=config_par)
    
    if args.post1d:
        update_key(config.pipe,'run_fibermap', False)
        update_key(config.pipe,'run_reduce2d', False)
        update_key(config.pipe,'run_combine', False)
        update_key(config.pipe,'run_fieldlist', False)
        update_key(config.pipe,'run_Summarymerge', False)
        update_key(config.pipe,'run_spcalib', False)
        update_key(config.pipe,'run_healpix', False)
        update_key(config.pipe,'skip2d', True)

    if config.pipe['email.allemail']:
        update_key(config.pipe, 'sendemail', True)

    if config.pipe['customSettings.custom_single_mjd']:
        update_key(config.pipe,'run_fibermap', False)
        update_key(config.pipe,'run_reduce2d', False)
        update_key(config.pipe,'run_combine', True)
        update_key(config.pipe,'run_analyze', True)
        update_key(config.pipe,'run_XCSAO',  True)
        update_key(config.pipe,'run_fieldlist', False)
        update_key(config.pipe,'run_fieldmerge', True)
        update_key(config.pipe,'run_specFiles', True)
        update_key(config.pipe,'run_Summarymerge', False)
        update_key(config.pipe,'run_spcalib', False)
        update_key(config.pipe,'run_healpix', False)
        update_key(config.pipe,'end2end', True)

    if config.pipe['customSettings.custom_name']:
        update_key(config.pipe,'custom', True)
        update_key(config.pipe,'run_fieldlist', False)
        update_key(config.pipe,'run_fibermap', False)
        update_key(config.pipe,'run_reduce2d', False)
        update_key(config.pipe,'run_healpix', False)
        update_key(config.pipe,'run_spcalib', False)
        update_key(config.pipe,'skip2d', True)
    else:
        update_key(config.pipe,'custom', False)
        if args.epoch:
            update_key(config.pipe,'run_reduce2d', False)
            update_key(config.pipe,'run_healpix', False)
            update_key(config.pipe,'run_fibermap', False)

    if args.sdssv:
        update_key(config.pipe,'MWM_fluxer', True)
        update_key(config.pipe,'noreject', True)
        update_key(config.pipe,'run_Summarymerge', False)
        update_key(config.pipe,'map3d', 'merge3d')

    update_key(config.pipe, 'obs',list(set(config.pipe['fmjdselect.obs'])) )

    if config.pipe['fmjdselect.epoch'] is None:
        update_key(config.pipe, 'epoch', False)

    if not config.pipe['fmjdselect.epoch']:
        update_key(config.pipe,'started', False)
        update_key(config.pipe,'abandoned', False)

    fill_none_with_false(config.pipe)
    fill_none_with_false(config.queue)



    if args.show_config:
        show_config()
        return

    queue = uubatchpbs()


