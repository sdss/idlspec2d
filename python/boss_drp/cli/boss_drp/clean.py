#!/usr/bin/env python3
from boss_drp.utils.cleanup_bkups import cleanup_bkups
from boss_drp.utils.clean import clean_fmjd
from boss_drp.field import field_to_string

from tqdm import tqdm

import click

@click.group(name='clean')
def clean():
    """Tools to clean the BOSS DRP outputs"""
    pass

@clean.command(name='backups')
@click.option("--topdir", envvar='BOSS_SPECTRO_REDUX', help="Boss Spectro Redux base directory")
@click.option("--run2d", envvar='RUN2D', help="Run2d")
@click.option("--epoch", is_flag=True, help="run for the epoch coadds")
@click.option("--custom", default=None, help="Name of custom Coadd")
@click.option("--backups", default=3, type=int, help="Number of backups to keep")
def clean_bkup(topdir, run2d, epoch, custom , backups):
    """ Clean the Summary Table File Backups"""
    if backups == 0:
        backups = None
    if backups is None:
        return
    
    cleanup_bkups(topdir, run2d, backups = backups, epoch=epoch,
                  custom = custom)

    
@clean.command(name='run')
@click.option("--clean_type", "--clean", required=True, 
              type=click.Choice(
                  ['all','spec2d','comb','spec1d','post','merge','reformat','spcalib'],
                  case_sensitive=False), help="Pipeline Step to start the cleaning")
@click.option("--topdir", envvar='BOSS_SPECTRO_REDUX',
              help="Optional override value for the environment variable $BOSS_SPECTRO_REDUX")
@click.option("--run2d",envvar='RUN2D', help="Optional override value for the environment variable $RUN2D")
@click.option("--run1d",envvar='RUN1D', help="Optional override value for the environment variable $RUN1D")
@click.option("--epoch",is_flag=True, help="Clean up epoch run")
@click.option("--reset",is_flag=True, help="if clean_type == all, then remove plans and redux")
@click.option("--remove_redux", is_flag=True, help="if clean_type == all, then remove redux")
@click.option("--dry",is_flag=True, help="Print Files to be removed rather then remove")
@click.option("--verbose",is_flag=True, help="Print Files paths (with wildcards) to be removed")
@click.option("--field", "-f", default='*', type=str, help="Run for a single Field")
@click.option("--mjd", "-m", default='?????', type=str, help="Run for a single MJD")
@click.option("--fmjd", multiple=True,help="List of Field-MJDs to clean")
def clean_run(clean_type,topdir,run2d,run1d,epoch,reset,remove_redux,dry,verbose,field,mjd,fmjd):
    """Clean pipeline prodcuts fro a given field, mjd, or field-mjd"""
    field = field_to_string(field)

    if not fmjd:
        fmjd = [f'{field}-{mjd}']
    for fmjd in tqdm(fmjd, leave=False, disable=not verbose):
        clean_fmjd(topdir, run2d, run1d, fmjd.split('-')[0],fmjd.split('-')[1],
                   epoch=epoch, dry=dry, clean_type=clean_type, reset=reset,
                   remove_redux = remove_redux, verbose = verbose )
