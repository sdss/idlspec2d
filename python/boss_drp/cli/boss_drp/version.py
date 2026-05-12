#!/usr/bin/env python3
try:
    from boss_drp import __version__
except ImportError:
    from sdsstools import get_package_version
    __version__ = get_package_version(__file__, 'boss_drp') or 'dev'

from boss_drp.utils.argparse_help import full_help_callback

import click

@click.command(name='version',context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.option(
    "--fullhelp",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=full_help_callback,
    help="Show full help including all subcommands")
def version():
    """Prints the IDLspec2D BOSS_DRP version"""
    click.echo(__version__)


if __name__ == "__main__":
    version()
