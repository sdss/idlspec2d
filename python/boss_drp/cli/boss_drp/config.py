from boss_drp.Config import prompt_edit, config as cfg_
from boss_drp import idlspec2d_dir
import click
from pathlib import Path
from ruamel.yaml import YAML
from sdsstools.configuration import DEFAULT_PATHS, get_config
import sys, os

@click.command(name='config')
@click.option('--load_file', default=None, help='Config File to Load as the base')
@click.option('--load_name', default='boss_drp', help='Config Name to Load as the base (eg. boss_drp,queue )')
@click.option('--queue', is_flag = True,  help='Set this flag if you are loading a queue config')
@click.option('--show', is_flag=True, default=False, help='Print the Config')
@click.option('--save', default=None, help='Location to Save the config to (required if you are editing)')
@click.option('--save_to_envvar', '-s',is_flag=True, help='Save to the ENVVAR (BOSS_DRP_PIPE_CONFIG_PATH or BOSS_DRP_PIPE_CONFIG_PATH) if set')
@click.option('--edit', is_flag=True, default=False, help='Edit the Config')
def config(load_file, load_name, queue, show, save,save_to_envvar, edit):
    """BOSS DRP Config Commands"""

    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.indent(mapping=2, sequence=4, offset=2)
    yaml.width = 100

    is_queue_file = load_file is not None and 'queue' in Path(load_file).name


    if queue or is_queue_file or load_name == 'queue':
        config_envvar = 'BOSS_DRP_QUEUE_CONFIG_PATH'
        if load_name != 'queue':
            click.echo(click.style("Warning: 'load_name' must be 'queue' for queue option/file.... setting", fg="yellow"))
            load_name = 'queue'
        load_file = os.path.join(idlspec2d_dir,'python','boss_drp','etc','queue.yml')
    else:
        config_envvar = 'BOSS_DRP_PIPE_CONFIG_PATH'


    if load_file is None:
        cfg_.load()
        load_name = cfg_.pipe._CONFIG_FILE

    else:
        try:
            _pipe = get_config(
                        load_name,
                        allow_user=True,
                        config_file=load_file,
                        config_envvar=config_envvar
                    )
        except:
            _pipe = get_config(
                        load_name,
                        allow_user=True,
                        config_file=load_file,
                        config_envvar=None
                    )
        load_name = _pipe._CONFIG_FILE

    
    with open(load_name, 'r') as f:
        click.echo(f'loading {load_name}')
        cfg = yaml.load(f)

    if save_to_envvar:
        if not save:
            save = os.getenv(config_envvar, None)
            if save:
                click.echo(f'Setting save to ${config_envvar}')

    if edit:
        if not save:
            raise click.ClickException("Save must be set to edit")
        cfg = prompt_edit(cfg)

    if show:
        yaml.dump(cfg, sys.stdout)

    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        click.echo(f'Saving to {save}')
        with open(save, "w") as f:
           yaml.dump(cfg, f)