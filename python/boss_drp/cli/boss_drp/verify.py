import click
import os
from boss_drp.utils.argparse_help import AttrDict, multi_str2bool, multi_str2none, full_help_callback, OrderedGroup

@click.command('verify', context_settings={"help_option_names": ['-h','--help'], "max_content_width": 150}) 
@click.option(
    "--fullhelp",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=full_help_callback,
    help="Show full help including all subcommands"
)
@click.argument("step", type=click.Choice(['setup','readfibermap','reduce','coadd','spec1d',
                                           'run_PyXCSAO','fieldlist','fieldmerge','spSpec_reformat',
                                           "spcalib_qa","all"],case_sensitive=False))
@click.option("--run2d", help='RUN2D version')
@click.option("--run1d", help='RUN1D version')
@click.option("--obs", help='Observatory of Coadd')
@click.option("--field", help='Field (or obs for custom coadds) of Coadd')
@click.option("--mjd", help='MJD of Coadd')
@click.option("--mjd1d", help='MJD of 1D RUN for Custom Coadds')
@click.option("--epoch", is_flag=True, help='Build for Epoch Coadds')
@click.option("--custom", help='Name of custom coadd')


@click.option("--skip", is_flag=True, help='Mark this step as complete')
@click.option("--verbose", is_flag=True, help='Print Status')
def verify(step,run2d, run1d, obs, field, mjd, mjd1d, epoch, custom, skip, verbose):
    """Verify BOSS DRP Steps"""
    from boss_drp.utils.verify_step import verify, setup_json

    if run2d is None:
        run2d = os.getenv('RUN2D')
    
    if step.lower() == 'setup':
        setup_json(run2d, field, mjd, mjd1d=mjd1d, epoch=epoch, 
                custom=custom, run1d=run1d, obs=obs)

    else:
        verify(step, run2d, field, mjd, epoch=epoch, custom=custom, mjd1d=mjd1d, 
               skip=skip, verbose=verbose)


