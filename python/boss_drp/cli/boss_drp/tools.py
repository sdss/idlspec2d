#!/usr/bin/env python3
from boss_drp.sos.sdR_hdrfix import fixhdr, getLastMJD
from boss_drp.prep import flag_manual_cal
from boss_drp.Flatlib.opfiber import refine_opfiber, build_trace_guess
from boss_drp.utils.argparse_help import AttrDict, full_help_callback, OrderedGroup
import click
from os import getenv
import platform

try:
    from termcolor import colored
except:
    def colored(text, color):
        return text

@click.group(name='tools', cls=OrderedGroup, context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.option(
    "--fullhelp",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=full_help_callback,
    help="Show full help including all subcommands"
)
def tools():
    """Miscellaneous BOSS DRP Tools"""
    pass



use_sos_style = ('sdss5' in platform.node()) or (getenv('IDLSPEC2D_SOS') is not None)
obs_env = (getenv('OBSERVATORY') or '').lower()

def get_cam_and_carts():
    if ('sdss5' in platform.node()) or (getenv('IDLSPEC2D_SOS') is not None):
        if getenv('OBSERVATORY', '').lower() == 'apo':
            return ['b1','r1','??'], ['FPS-N']
        else:
            return ['b2','r2','??'], ['FPS-S']
    else:
        return ['b1','b2','r1','r2','??'], ['FPS-S', 'FPS-N']

cam_choices, cart_choices = get_cam_and_carts()


def hdrfix_cmds():
    def decorate(fn):
        fn = click.argument('expid', metavar='EXPID')(fn) 
        fn = click.option('--mjd', '-m', default=None, help='MJD of file (default: latest)')(fn)

        req = True if ((obs_env == '') or (obs_env is None)) else False
        hidden = True if use_sos_style else False
        fn = click.option('--obs', type=click.Choice(['APO', 'LCO'], case_sensitive=False),
                          envvar='OBSERVATORY', required=req, help='Observatory', hidden = hidden, )(fn)

        fn = click.option('--clobber', is_flag=True, help='Clobber sdHdrFix file')(fn)
        fn = click.option('--cameras', type=click.Choice(cam_choices, case_sensitive=False), 
                          default='??', help='Cameras for hdr update')(fn)
        fn = click.option('--no-update', '-u', is_flag=True, help='Skip updating SOS logs')(fn)
        fn = click.option('--nogit', is_flag=True, help='Skip automatic git add')(fn)

        # Quality flags
        fn = click.option('--bad', '-b', is_flag=True, help='Flag as bad')(fn)
        fn = click.option('--test', '-t', is_flag=True, help='Flag as test')(fn)

        # Lamp options
        fn = click.option('--FF', nargs=4, type=click.Choice(['0','1']),help='Flat Field Lamp')(fn)
        fn = click.option('--FFS', nargs=8, type=click.Choice(['0','1']),help='Flat Field Screen')(fn)
        fn = click.option('--NE', nargs=4, type=click.Choice(['0','1']),help='Ne arc lamp')(fn)
        if use_sos_style:
            if obs_env == 'apo':
                fn = click.option('--HGCD', nargs=4, type=click.Choice(['0', '1']), help='HgCd arc Lamp')(fn)
            else:
                fn = click.option('--HEAR', nargs=4, type=click.Choice(['0', '1']), help='HeAr arc Lamp')(fn)
        else:
                fn = click.option('--HGCD', nargs=4, type=click.Choice(['0', '1']), help='HeCd arc Lamp')(fn)
                fn = click.option('--HEAR', nargs=4, type=click.Choice(['0', '1']), help='HeAr arc Lamp')(fn)


        fn = click.option('--arc', is_flag=True, default=False, 
                          help='short cut to set all relevant arc lamps to 1 1 1 1')(fn)
        fn = click.option('--flat', is_flag=True, default=False, 
                          help = 'short cut to set FF = 1 1 1 1 & FFS =  1 1 1 1 1 1 1 1 ')(fn)
        fn = click.option('--hartmann', help='Hartmann Door Status',
                          type=click.Choice(['Out', 'Right', 'Left', 'Closed'], case_sensitive=False))(fn)

        # Common
        fn = click.option('--quality', help='Set Quality flat of exposures',
                          type=click.Choice(['excellent', 'test', 'bad'], case_sensitive=False))(fn)

        # Specialized
        fn = click.option('--flavor',
              type=click.Choice(['bias','dark','flat','arc','science','smear']),help='Type/Flavor of exposure')(fn)
        fn = click.option('--exptime', type=float, help='Exposure length (s)')(fn)
        fn = click.option('--tai-beg', 'TAI_BEG', type=float, help='Starting time (tai) of exposure')(fn)
        fn = click.option('--cartid', type=click.Choice(cart_choices), help='Cartridge Mounted')(fn)
        fn = click.option('--fieldid', type=int, help="FieldID")(fn)
        fn = click.option('--confid', type=int, help="ConfigureID")(fn)
        fn = click.option('--designid', type=int, help='DesignID')(fn)

        # Manual key
        fn = click.option('--key', '-k', help = 'header keyword to update (required if value is set)')(fn)
        fn = click.option('--value', '-v', help = 'updated header keyword value (required if key is set)')(fn)
        return fn
    return decorate

@tools.command(name='sdR_hdrfix',context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@hdrfix_cmds()
@click.pass_context
def run_sdR_hdrfix(ctx, **kwrds):
    """Create the files used by the pipeline to fix the header meta data of the BOSS exposures"""
    args = AttrDict(ctx.params)

    if args.mjd is None:
        args.mjd = str(getLastMJD(silent=True))

    if args.obs is None:
        args.obs = getenv('OBSERVATORY')

    if not args.obs:
        raise click.ClickException('OBSERVATORY is required')

    args.obs = args.obs.lower()

    # shortcuts
    if args.arc:
        if args.obs.lower() == 'apo':
            args.ne = ('1','1','1','1')
            args.hgcd = ('1','1','1','1')
        else:
            args.ne = ('1','1','1','1')
            args.hear = ('1','1','1','1')

    if args.flat:
        args.ff = ('1','1','1','1')
        args.ffs = ('1','1','1','1','1','1','1','1')

    # build updates dict
    updates = {}

    if args.bad:
        updates = {'quality': 'bad'}
    elif args.test:
        updates = {'quality': 'test'}
    else:
        skip = {'expid','mjd','obs','clobber','bad','test','key','value',
                'rerun','camera','observer','nogit','no_update','arc','flat'}

        for k, v in args.items():
            if k in skip:
                continue
            if v is None:
                continue
            if k.upper() in {'DATE_OBS', 'TAI_BEG'}:
                k = k.replace('_','-')
            if (isinstance(v, tuple) or isinstance(v, list)):
                updates[k] = ' '.join(v)
            else:
                updates[k] = str(v)

        if args.key and args.value:
            if args.key.upper() in ['EXPOSURE', 'MJD']:
                click.echo('ERROR: Invalid Keys. Exiting')
                return
            updates[args.key] = args.value

    # confirmation loop
    for k in list(updates.keys()):
        if k.lower() in {'fieldid','confid','designid','flavor','exptime','tai-beg','cartid'}:
            confirm = click.prompt(f'Do you really want to edit {k}? (y/[n])', default='n')
            if confirm.lower() != 'y':
                updates.pop(k)
                click.echo(f'Skipping {k}')

    if not updates:
        click.echo('No updates specified.')
        return

    fixhdr(args.expid, updates,
           mjd=args.mjd,
           obs=args.obs,
           clobber=args.clobber,
           cameras=args.cameras,
           update=(not args.no_update),
           nogit=args.nogit)



@tools.command(name='flag_manual_cal', context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.option("-o", "--observatory", "--obs", required=True, help="Observatory",
              type=click.Choice(["apo", "lco"], case_sensitive=False))
@click.option("-m", "--mjd",required=True, type=int, help="MJD")
@click.option("-f", "--field", required=True, type=str, help="FieldID")
@click.option("-e", "--expid", required=False, type=int, default=None,
              help="Exposure ID to manually set the calibration frame exposure ID")
@click.option("-t", "--type", "cal_type", required=True, help="Calibration Type",
              type=click.Choice(["arc", "flat"], case_sensitive=False))
@click.option("--nogit", is_flag=True, help="Skip automatic git add")
def run_flag_manual_cal(observatory, mjd, field, expid, cal_type, nogit):
    """Build spManCal.par file to flag manual alternative calibration frames for spPlan"""
    flag_manual_cal(type=cal_type, field=field, mjd=mjd,
                    obs=observatory.lower(), expid=expid, nogit=nogit)



@tools.group(name='opFiber')
@click.option(
    "--fullhelp",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=full_help_callback,
    help="Show full help including all subcommands"
)
def opFiber():
    """Prints updated/refined values for opfibers using spFlat traces"""
    pass

@opFiber.command(name = 'Refine')
@click.argument("fitsfile",type=click.Path(exists=True, dir_okay=False, readable=True))
@click.option("-p", "--precision", type=int, default=3,show_default=True,
              help="Precision of the reported fiberspacing and bundle gaps.")
def refine(fitsfile, precision):
    """Refine the opFiberFPS parameters using a spFlat."""
    refine_opfiber(fitsfile, precision=precision)


@opFiber.command(name='Guess', short_help='Guess at the opFiberFPS parameters using a sdProc-XX-XXXXXXXX.fits file',
                 context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.argument("sdProcFile", type=click.Path(exists=True, dir_okay=False, readable=True))
@click.option("-b", "--bundlefibers", multiple=True, type=int,
              help=("List of number of fibers per bundle. Use multiple times, for example: -b 2 -b 4 -b 4"))
@click.option("-m", "--mjd", type=int, help="MJD of new updated OpFiber Fiberparameter entry.")
@click.option("-f", "--plot",is_flag=True,default=False,  
              help="Whether to plot the flux slice and detected peaks.")
@click.option("--min_peak_sep",type=float, default=6, show_default=True, 
              help="Minimum separation between detected peaks.")
@click.option("--min_peak_height", type=float, default=5000, show_default=True, 
              help="Minimum flux level to be detected as a peak.")
@click.option("-p","--precision", type=int,default=3,show_default=True, 
              help="Precision of the reported fiberspacing and bundle gaps.")
def guess(sdProcFile, bundlefibers, mjd, plot, min_peak_sep, min_peak_height, precision):
    """Takes a processed image frame produced using the /sawraw flag in sdssproc (or indirecly via spreduce2d)
       as an input. It then uses either the bundlefiber list of number of fibers per bundle
       (supplied as input or via opFiberFPS) combined with the scipy peak finding algarithm
       to create a first guess of the peak fiber and bundle gaps. It uses the median
       flux of the 11 central pixel (along the dispersion axis) to build the flux array"""
    build_trace_guess( sdProcFile, bundlefibers=list(bundlefibers) if bundlefibers else None, mjd=mjd,
                      plot=plot, min_peak_sep=min_peak_sep, min_peak_height=min_peak_height, precision=precision )



if __name__ == '__main__':
    tools()