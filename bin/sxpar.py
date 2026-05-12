#!/usr/bin/env python3

import click
from boss_drp.utils.sxpar import sxpar
from boss_drp.utils.argparse_help import full_help_callback


@click.command(context_settings=dict(help_option_names=['-h', '--help'],max_content_width= 150)) 
@click.option(
    "--fullhelp",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=full_help_callback,
    hidden = True,
    help="Show full help including all subcommands"
)
@click.argument("fitsfile", type=click.Path(exists=True, dir_okay=False, readable=True))
@click.argument("keyword")
@click.option("-v", "--verbose", is_flag=True, help="verbose")
def cli(fitsfile, keyword, verbose):
    """Simply parse a fits header."""
    output = sxpar(fitsfile, keyword, verbose)
    for line in output:
        click.echo(line)


if __name__ == "__main__":
    cli()