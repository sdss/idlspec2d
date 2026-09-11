import click
from .boss_drp.plan import plan
from .boss_drp.config import config
from .boss_drp.tools import tools
from .boss_drp.log import log
from .boss_drp.clean import clean
from .boss_drp.daily import daily
from .boss_drp.batch import batch
from .boss_drp.run import run
from .boss_drp.version import version
from .boss_drp.verify import verify
from boss_drp.utils.argparse_help import full_help_callback, OrderedGroup

@click.group(cls=OrderedGroup, context_settings={"help_option_names": ["-h", "--help"],"max_content_width": 150})
@click.option(
    "--fullhelp",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=full_help_callback,
    help="Show full help including all subcommands"
)


def cli():
    pass

for cmd in (version, config, plan, daily, batch, run, log, clean, tools, verify):
    cli.add_command(cmd)


if __name__ == "__main__":
    cli()