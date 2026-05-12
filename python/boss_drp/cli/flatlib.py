from boss_drp.Flatlib import (
    analysis as flat_analysis,
    build as flat_build,
    plot as flat_plot,
    reduce as flat_reduce,
    read_fiberAssignments,
)
from boss_drp import QA_DIR
from boss_drp.utils import jdate
from boss_drp.Config import config, fill_none_with_false

from boss_drp.utils.argparse_help import full_help_callback, OrderedGroup
import click
import os.path as ptt
import glob
from os import getenv, environ, makedirs

# ---- Shared options (like parent parsers) ----
def common_options(f):
    f = click.option('--dir', '-d', 'dir_', default=None, help='Flat Library Directory', metavar='FLATLIB_DIR')(f)
    f = click.option('--run2d', default=None, help='Override $RUN2D', metavar="RUN2D")(f)
    f = click.option('--lco/--no-lco', default=None, flag_value=True, help='Run for LCO data')(f)
    return f

def queue_options(f):
    f = click.option('--link/--no-link', is_flag=True, default=False, help='Link Pre-existing spFlat Files')(f)
    f = click.option('--deep/--no-deep', is_flag=True, default=False, help='Check Pre-existing plans for completion')(f)
    f = click.option('--queue_config', '--qc', default='flatlib', help='Queue Config name', metavar='QUEUE_CONFIG_NAME')(f)
    f = click.option('--queue_config_file', '--qcf', default=None, help='Queue Config File Path', metavar='QUEUE_CONFIG_PATH')(f)
    f = click.option('--nodes', type=int, default=None, help='Number of nodes to use', metavar='NNODES')(f)
    f = click.option('--submit/--no-submit', is_flag=True, default=None, help='Submit the job to the queue')(f)
    f = click.option('--run', is_flag=True, default=False, help='Just link (if set), but do not run new spFlat files')(f)
    f = click.option('--link_all', '-a', is_flag=True, default=False, help='Link all spFlat files regardless of spPlanTrace file')(f)
    f = click.option('--link_traceflat', '-c',is_flag=True, default=False, help='Link all spTraceFlat files')(f)
    return f

def timeSeries_options(f):
    f = click.option('--mjd', '-m', multiple=True, help='List of mjds to plot alone', metavar='MJD')(f)
    f = click.option('--mjdstart', type=int, default=None, help='MJD to start reduction', metavar='MJDSTART')(f)
    f = click.option('--TraceIDs', '-t',is_flag=True,  default=False,
                     help='Label with Trace FiberIDs rather then slit FiberIDs')(f)
    return f

def resolve_common_path(dir_value, run2d):
    if run2d is not None:
        environ["RUN2D"] = run2d

    if dir_value is None:
        dir_value = ptt.join(QA_DIR, "BOSSFLATLIB", getenv("RUN2D"))

    makedirs(dir_value, exist_ok=True)
    return dir_value


def load_queue_if_needed(command, queue_config, queue_config_file, nodes, nosubmit):
    if command in ("end2end", "reduce"):
        config.load(queue_config, queue_config_file=queue_config_file)
        config_par = {
            "nodes": nodes,
            "no_submit": nosubmit,
        }
        for name, val in config_par.items():
            config.queue.set(name, val)

        fill_none_with_false(config.queue)



# ---- Main CLI group ----

@click.group(cls = OrderedGroup, context_settings={"help_option_names": ['-h','--help'], "max_content_width": 150}) 
@click.option(
    "--fullhelp",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=full_help_callback,
    help="Show full help including all subcommands"
)
@click.pass_context
def cli(ctx, **kwrds):
    """Build and analyze a library of flats to check for Fiber throughput Issues"""
    ctx.ensure_object(dict)


# ---- Commands ----

@cli.command(name='reduce', context_settings={"help_option_names": ["-h", "--help"],"max_content_width": 150})
@common_options
@queue_options
@click.option('--mjd', '-m', multiple=True, help='MJDs to Run', metavar='MJD')
@click.option('--fps', is_flag=True, default=False, help='Catch up FPS')
@click.option('--plates', is_flag=True, default=False, help='Catch up Plates')
@click.option('--legacy', is_flag=True, default=False, help='Catch up Legacy')
@click.pass_context
def reduce_cmd(ctx, dir_, run2d, lco, mjd, fps, plates, legacy, link, deep, queue_config, queue_config_file,
               nodes, submit, run, link_all, link_traceflat):
    """Reduce/link the spFlats"""
    load_queue_if_needed("reduce", queue_config, queue_config_file, nodes, submit)

    flat_reduce(
        resolve_common_path(dir_, run2d),
        list(mjd) if mjd else None,
        link=link,
        deep=deep,
        lco=lco,
        plates=plates,
        legacy=legacy,
        fps=fps,
        link_all=link_all,
        nosubmit=submit,
        nodes=nodes,
        no_run=run,
        link_traceflat=link_traceflat,
    )


@cli.command(name='build', context_settings={"help_option_names": ["-h", "--help"],"max_content_width": 150})
@common_options
@click.pass_context
def build(ctx, dir_, run2d, lco):
    """Build the flat library fits file"""
    obs = ["lco"] if lco else ["apo"]
    flat_build(resolve_common_path(dir_, run2d), obs)


@cli.command(name='plot', context_settings={"help_option_names": ["-h", "--help"],"max_content_width": 150})
@common_options
@click.option('--save', '-s', default=None, help='Save Directory', metavar='SAVEDIR')
@click.option('--mjd', '-m', multiple=True, help='List of mjds to plot', metavar='MJD')
@click.option('--flats', '-f', multiple=True, help='List of reduced flats to plot', metavar='FLATLIST')
@click.pass_context
def plot(ctx, dir_, run2d, lco, save, mjd, flats):
    """Plot Raw and Reduced Flat"""
    obs = ["lco"] if lco else ["apo"]
    dir_ = resolve_common_path(dir_, run2d)
    if mjd:
        mjd_list = list(mjd)
    else:
        mjd_list = []
        for ob in obs:
            mjd_list.extend(
                [ptt.basename(x) for x in glob.glob(ptt.join(dir_, "calibs", ob.lower(), "*"))]
            )
        mjd_list = list(set(mjd_list))

    assigns = {}
    for ob in obs:
        assigns[ob] = read_fiberAssignments(
            ptt.join(dir_, "fiberAssignments", ob, "fiberAssignments.csv")
        )

        if flats:
            flat_list = []
            for flat in flats:
                flat_list.extend(glob.glob(ptt.join(dir_, "calibs", ob.lower(), "*", flat)))
        else:
            flat_list = []
            for one_mjd in mjd_list:
                flat_list.extend(
                    glob.glob(ptt.join(dir_, "calibs", ob.lower(), one_mjd, "spFlat*.fits*"))
                )

        save_dir = save if save else ptt.join(dir_, "plots", ob)

        for flat in flat_list:
            flat_plot(dir_, flat, savedir=save_dir, assigns=assigns)


@cli.command(name="analyze", context_settings={"help_option_names": ["-h", "--help"],"max_content_width": 150})
@common_options
@click.option('--mjd', '-m', multiple=True, help='List of mjds to plot alone', metavar='MJD')
@click.option('--plot', default=False, is_flag=True, help='Plot Flat')
@click.pass_context
def analyze_cmd(ctx, dir_, run2d, lco, mjd, plot):
    """Run Full analysis on Flat library"""
    obs = "lco" if lco else "apo"
    dir_ = resolve_common_path(dir_, run2d)
    parms = {"obs": obs, "run": "all", "noplot": plot}

    if mjd:
        for one_mjd in mjd:
            flat_analysis(dir_, getenv('RUN2D'), mjd=one_mjd, **parms)
    else:
        flat_analysis(dir_, getenv('RUN2D'), **parms)


@cli.command(name='lowfiber', context_settings={"help_option_names": ["-h", "--help"],"max_content_width": 150})
@common_options
@click.option('--mjd', '-m', multiple=True, help='List of mjds to plot alone', metavar='MJD')
@click.option('--threshold', '-t', type=float, default=0.8, help='Threshold to flag lowfibers', metavar='THRESHOLD')
@click.pass_context
def lowfiber(ctx, dir_,run2d, lco, mjd, threshold):
    """Check for Low fibers"""
    dir_ = resolve_common_path(dir_, run2d)
    obs = "lco" if lco else "apo"
    parms = {"obs": obs, "run": "lowfiber", "lowFiber": threshold}

    if mjd:
        for one_mjd in mjd:
            flat_analysis(dir_, getenv('RUN2D'), mjd=one_mjd, **parms)
    else:
        flat_analysis(dir_, getenv('RUN2D'), **parms)


@cli.command(name='csv', context_settings={"help_option_names": ["-h", "--help"],"max_content_width": 150})
@common_options
@click.pass_context
def csv(ctx, dir_, run2d, lco):
    """Export CSV only"""
    obs = "lco" if lco else "apo"
    dir_ = resolve_common_path(dir_, run2d)

    flat_analysis(dir_, getenv('RUN2D'), obs=obs, run="csv")


@cli.command(name="timeSeries", context_settings={"help_option_names": ["-h", "--help"],"max_content_width": 150})
@common_options
@timeSeries_options
@click.pass_context
def timeseries(ctx, dir_, run2d, lco, mjd, mjdstart, traceids):
    """Plot Throughout Time Series only"""
    obs = "lco" if lco else "apo"
    parms = {"obs": obs, "run": "timeSeries", "TraceIDs": traceids}
    dir_ = resolve_common_path(dir_,run2d)
    if mjd:
        for one_mjd in mjd:
            flat_analysis(dir_, getenv("RUN2D"), mjd=one_mjd, **parms)
    else:
        flat_analysis(dir_, getenv("RUN2D"), **parms)


@cli.command(name='end2end', context_settings={"help_option_names": ["-h", "--help"],"max_content_width": 150})
@common_options
@queue_options
@timeSeries_options
@click.pass_context
def end2end(ctx, dir_, run2d, lco, mjd, mjdstart, link, deep, queue_config, queue_config_file, nodes,
                submit, run, link_all, link_traceflat, traceids):
    """Run full pipeline and plot time series (FPS only)"""
    load_queue_if_needed("end2end", queue_config, queue_config_file, nodes, submit)

    dir_ = resolve_common_path(dir_,run2d)
    obs = "lco" if lco else "apo"

    mjd_list = list(mjd) if mjd else []
    if mjdstart is not None:
        if mjdstart < 0:
            mjd_list.append(jdate.astype(int) + mjdstart)

    flat_reduce(
        dir_, mjd_list,
        link=link, deep=deep, lco=lco,
        plates=False, legacy=False, fps=True, link_all=link_all,
        nosubmit=submit, nodes=nodes, no_run=run,
        link_traceflat=link_traceflat, mjdstart=mjdstart
    )

    flat_build(dir_, [obs])

    parms = {"obs": obs, "noplot": True, "run": "csv", "TraceIDs": traceids}
    flat_analysis(dir_, getenv("RUN2D"), **parms)

    parms["run"] = "timeSeries"
    flat_analysis(dir_, getenv("RUN2D"), **parms)



# ---- Entry point ----
if __name__ == "__main__":
    cli()

