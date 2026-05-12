from boss_drp.run.uurundaily import uurundaily
from boss_drp.utils.argparse_help import parse_num_list, AttrDict
from boss_drp.Config import config, update_key, show_config, show_config_opt, fill_none_with_false
from boss_drp.cli.boss_drp.cli2config import cli2config

import os
import click


@click.command(name = 'daily',context_settings={"help_option_names": ["-h", "--help"]}, short_help='Run Pipeline Planning to Post')
@click.option("--module", type=str, default=None, help="Module for daily run")

# Field-MJD Selection
@click.option("--apo", "obs", flag_value="apo", multiple=True, help="Run apo")
@click.option("--lco", "obs", flag_value="lco", multiple=True, help="Run lco")
@click.option("--mjd", type=int, multiple=True,
              help="Manually run for a single/list of mjd (does not update nextmjd.par)")
@click.option("--range_mjd","range_mjd", type=str,
              callback=lambda ctx, param, value: parse_num_list(value) if value is not None else None,
              help="Manually run for a range of mjds (does not update nextmjd.par)")
@click.option("--dither/--no-dither", "dither", is_flag=True, default=None,
              help="Skip Dither Engineering Fields")
@click.option("--epoch/--no-epoch", is_flag=True, default=None, help="Run Epoch Coadds")
@click.option("--increment_nextmjd","--nextmjd", "increment_nextmjd", is_flag=True, default=None, 
              help='Increment NextMJD file even when specifying an MJD')
@click.option("--daily_epoch" "--de", "daily_epoch", is_flag = True, default = None, 
              help="Run Epoch Coadds after daily coadds")

# Pipeline Steps
@click.option("--traceflat/--no-traceflat","run_spTrace", is_flag=True, default=None,
              help="Skip Building and using TraceFlats")
@click.option("--force-arc2trace/--no-force-arc2trace", "force_arc2trace", is_flag=True, default=None,
              help="Force Use of Arc2Trace")
@click.option("--fibermap/--no-fibermap", "run_fibermap", is_flag=True, default=None,
              help="Skip Pre-Run of readfibermap")
@click.option("--no-prep/--prep", "no_prep", is_flag=True, default=None,
              help="Skip building TraceFlats and spfibermaps before pipeline run")
@click.option("--skip-plan", "skip_plan", multiple = True, default = (), help="Skip the given plan",
              type=click.Choice(["pipe", "trace", "all"], case_sensitive=False))
@click.option("--skip-plan-all", is_flag=True, default=False, help="Skip all plans")
@click.option("--clobber", "clobber", multiple = True, default= (),
               type=click.Choice(["spPlans", "fibermap", "trace", "pipe", "all"], case_sensitive=False),
               help="Clobber uubatchpbs + a combo of spPlan, fibermap, and TraceFlat run")
@click.option("--clobber-all", "clobber_all", is_flag=True, help="Clobber everything")

@click.option("--healpix/--no-healpix","run_healpix", is_flag=True, default=None,
              help="Turn off copy to healpix")
@click.option("--summary/--no-summary", is_flag=True, default=None, help="Build Summary Files")

# Debug
@click.option("--saveraw/--no-saveraw", is_flag=True, default=None,help="save sdssproc outputs")
@click.option("--debug/--no-debug", is_flag=True, default=None,help="save extraction debug files")

# Shortcut Options
@click.option("--tagged", "tagged_run", is_flag=True, default=False,
              help="sets --merge3d --sc tagged_daily --no-dither --monitor --allemail")
@click.option("--daily", "daily_run", is_flag=True, default=False,
              help="sets --merge3d --sc fast_daily --monitor --allemail --no-healpix")
@click.option("--dev", "dev_run", is_flag=True, default=False,
              help="sets --merge3d --sc tagged_daily --no-dither --monitor --no-healpix")

# Pipeline Options
@click.option("--topdir","--top-dir","BOSS_SPECTRO_REDUX", type=str, default=None,
              help="Optional override value for the config") 
@click.option("--run1d", "RUN1D", type=str, default=None,
              help="Optional override value for config") 
@click.option("--run2d", "RUN2D", type=str, default=None,
              help="Optional override value for the config")
@click.option("--nodist/--dist", "nodist", is_flag=True, default=None,
              help="unsets/sets --nodist and reactivates/deactivates the flux distortion corrections")
@click.option("--bay15", "map3d", flag_value="bayestar15", default=None, help="Set map3d to bayestar15 model")
@click.option("--merge3d", "map3d", flag_value="merge3d", default=None, help="Set map3d to best 3d model")
@click.option("--batch/--no-batch", "batch_mjd", is_flag=True, default=None, help="run for multiple mjds in a single batch")
@click.option("--nodb/--db", "no_db", is_flag=True, default=None, help="skip Database operations")
@click.option("--monitor/--no-monitor", "pipe_monitor", is_flag=True, default=None, help="Monitors pipeline status")
@click.option("--pause", type=int, default=None, help="Pause time (s) in status updates")
@click.option("--allemail/--no-allemail",  is_flag=True,  default=None,
              help="Email intermediate log using all emails in $DAILY_DIR/etc/emails (defaults to first email only)")

# Configurations
@click.option("--pipe_config", "--pipe-config", "--pc", "config", default="boss_drp", help="Queue Config name")
@click.option("--pipe_config_file", "--pipe-config-file", "--pcf", "config_file", default=None, help="Queue Config File Path")
@click.option("--queue_config", "--queue-config", "--qc", "queue_config", default=None, help="Queue Config name")
@click.option("--queue_config_file", "--queue-config-file", "--qcf", "queue_config_file", default=None, help="Queue Config File Path")

# Queue Options
@click.option("--no-write/--write", "no_write", is_flag=True, default=None, help="skip writing and submitting job")
@click.option("--nosubmit/--submit", is_flag=True, default=None, help="Skip submitting uubatch job (ideal for allowing editting of plans)")
@click.option("--walltime", help="Wall time in hours", type=str, default=None)
@click.option("--mem_per_cpu", "--mem-per-cpu", help="Memory allocated per CPU", type=str, default=None)
@click.option("--nbundle", help="Number of jobs to bundle", type=int, default=None)
@show_config_opt
@click.pass_context
def daily(ctx, **kwrds):
    """Plan, run Spectro-2D and Spectro-1D reductions, and run post pipeline steps"""
    args = AttrDict(ctx.params)
    # argparse MultiBoolAction equivalent
    
    if args['no_prep']:
        args.run_fibermap = False
        args.run_spTrace = False

    # Normalize tuple outputs into lists for easier downstream use
    args.obs = list(args.obs) if args.obs else None
    args.mjd = list(args.mjd) if args.mjd else None
    args.clobber = list(args.clobber or ())

    if args.BOSS_SPECTRO_REDUX is not None:
        os.environ['BOSS_SPECTRO_REDUX'] = args.BOSS_SPECTRO_REDUX
    if args.RUN2D is not None:
        os.environ['RUN2D'] = args.RUN2D
    if args.RUN1D is not None:
        os.environ['RUN1D'] = args.RUN1D

    if (args.tagged_run) or (args.daily_run) or (args.dev_run):
        args.map3d = 'merge3d'
        args.pipe_monitor = True
    
    if args.tagged_run:
        args.dither = False
        args.allemail = True
        if args.queue_config is None:
            args.queue_config = 'tagged_daily'
    elif args.daily_run:
        args.allemail = True
        args.run_healpix = False
        if args.queue_config is None:
            args.queue_config = 'fast_daily'
    elif args.dev_run:
        args.dither = False
        args.run_healpix = False
        if args.queue_config is None:
            args.queue_config = 'tagged_daily'
    elif args.queue_config is None:
        args.queue_config = 'tagged_daily'

    config_par = {'mem_per_cpu':args.mem_per_cpu,
                  'wall':args.walltime,
                  'no_submit': args.nosubmit,
                  'no_write': args.no_write,
                  'nbundle':args.nbundle}
    
    exclude_args = ['no_write','nosubmit','walltime','mem_per_cpu','nbundle',
                    'pipe_config','pipe_config_file','queue_config','queue_config_file',
                    'dev_run', 'skip_plan_all', 'clobber_all', 'tagged_run', 'daily_run', 
                    'show_config','clobber']
    cli2config(args, config_par =config_par, exclude=exclude_args, set_gen=False)


    if args.range_mjd is not None:
        if config.pipe['fmjdselect.mjd'] is not None:
            if not isinstance(config.pipe['fmjdselect.mjd'], list):
                update_key(config.pipe, 'mjd', list(config.pipe['fmjdselect.mjd']))
            config.pipe['fmjdselect.mjd'] + update_key(config.pipe, 'mjd', list(config.pipe['fmjdselect.mjd']))
 

        else:
            update_key(config.pipe, 'mjd', list( args.range_mjd))

    if args.clobber_all:
        args.clobber = ("all",)

    args.clobber = list(args.clobber or ())
    for c in args.clobber:
        if c.lower() == 'fibermap':
            update_key(config.pipe, 'clobber_pipe', True)
            update_key(config.pipe, 'clobber_fibermap', True)
        elif c.lower() == 'trace':
            update_key(config.pipe, 'clobber_pipe', True)
            update_key(config.pipe, 'clobber_spTrace', True)
        elif c.lower() == 'pipe':
            update_key(config.pipe, 'clobber_pipe', True)
        elif c.lower() == 'spplans':
            update_key(config.pipe, 'clobber_plan', True)
        elif c.lower() == 'all':
            for key in config.pipe.get('Clobber').keys():
                update_key(config.pipe, key, True)

    args.skip_plan = list(args.skip_plan or ())
    if args.skip_plan_all:
        args.skip_plan = ["all"]

    for sp in args.skip_plan:
        if sp.lower() == 'pipe':
            update_key(config.pipe, 'run_plan', False)
        elif sp.lower() == 'trace':
            update_key(config.pipe, 'run_spTrace_plan', False)
        elif sp.lower() == 'all':
            update_key(config.pipe, 'run_plan', False)
            update_key(config.pipe, 'run_spTrace_plan', False)

    if args.summary:
        update_key(config.pipe, 'run_Summarymerge', True)

    if args.module is None:
        args.module = os.getenv('MODULE', default=None)
        if args.module is None:
            args.module = os.getenv('RUN2D', default=None)
            if args.module is None:
                args.module = 'bhm/master'
            elif args.dev_run:
                args.module = f'work/{args.module}'
            else:
                args.module = f'bhm/{args.module}'
    update_key(config.pipe, 'module', args.module)
                
    if not args.obs:
        update_key(config.pipe, 'obs',['apo','lco'] )
        

    if config.pipe['general.V_TARG'] is None:
        update_key(config.pipe, 'V_TARG', '*')

    if config.pipe['reduce.map3d'] is None:
        update_key(config.pipe, 'map3d', 'merge3d')
    if config.pipe['monitor.pause'] is None:
        update_key(config.pipe, 'pause', 15*60)

    if config.pipe['fmjdselect.epoch'] is None:
        update_key(config.pipe, 'epoch', False)

    if not config.pipe['fmjdselect.epoch']:
        update_key(config.pipe,'started', False)
        update_key(config.pipe,'abandoned', False)

    fill_none_with_false(config.pipe)
    fill_none_with_false(config.queue)

    if args.show_config:
        show_config()
    else:
        uurundaily()

if __name__ == "__main__":
    daily()