
import argparse
import click
import sys
import re

class FullHelpAction(argparse.Action):
    def __init__(self, option_strings, dest=argparse.SUPPRESS, default=argparse.SUPPRESS, help=None):
        super(FullHelpAction, self).__init__(option_strings=option_strings, dest=dest, default=default, nargs=0, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        print_full_help(parser)
        sys.exit()
        
def print_full_help(parser):
    """Print full help including all subparsers."""
    parser.print_help()
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            for choice, subparser in action.choices.items():
                print(f"\nSubparser '{choice}' help:")
                subparser_help = subparser.format_help()
                print(''.join(' '*4 + line for line in subparser_help.splitlines(True)))
    sys.exit()

# class MultiBoolAction(argparse.Action):
#     def __init__(self, option_strings, dest, dests=None, expose_self=True, **kwargs):
#         self.dests = dests or []
#         self.expose_self = expose_self
#         #kwargs.setdefault("nargs", 0)
#         kwargs.setdefault("default", argparse.SUPPRESS)
#         super().__init__(option_strings, dest, nargs=0, **kwargs)
#     def __call__(self, parser, namespace, values, option_string=None):
#         for d in self.dests:
#             setattr(namespace, d, self.const)

#         # remove the temporary attribute
#         if not self.expose_self:
#             if hasattr(namespace, self.dest):
#                 delattr(namespace, self.dest)


class MultiBoolAction(argparse.Action):
    def __init__(self, option_strings, dest, dests=None, expose_self=False, **kwargs):
        self.dests = dests or []
        self.expose_self = expose_self
        kwargs.setdefault("nargs", 0)
        kwargs.setdefault("default", argparse.SUPPRESS)
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        for d in self.dests:
            setattr(namespace, d, False)

        if self.expose_self:
            setattr(namespace, self.dest, True)



def print_full_help_click(cmd: click.Command, ctx: click.Context, indent: int = 0) -> None:
    click.echo(cmd.get_help(ctx))

    if isinstance(cmd, click.Group):
        for name in cmd.list_commands(ctx):
            subcmd = cmd.get_command(ctx, name)
            if subcmd is None:
                continue

            click.echo(f"\n{' ' * indent}Subcommand '{name}' help:")
            subctx = click.Context(subcmd, info_name=name, parent=ctx)
            subhelp = subcmd.get_help(subctx)
            click.echo("".join(" " * (indent + 4) + line for line in subhelp.splitlines(True)))

            if isinstance(subcmd, click.Group):
                print_full_help_click(subcmd, subctx, indent + 4)


def full_help_callback(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    print_full_help_click(ctx.command, ctx)
    ctx.exit()

def multi_bool_option(*param_decls, dests=None, expose_self=False, **kwargs):
    return click.option(
        *param_decls,
        is_flag=True,
        callback=print_full_help_click(dests, expose_self),
        **kwargs
    )


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        for key, value in list(self.items()):
            if isinstance(value, tuple):
                value = list(value)
                if len(value) == 0:
                    value = None
                self[key] = value

    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


    
def parse_num_list(ctx, param, value):
    """
    Replace this with your real parseNumList logic.
    Return a list of ints or whatever your app expects.
    """
    if value is None:
        return [None]

    m = re.match(r'(\d+)(?:-(\d+))?$', value)
    # ^ (or use .split('-'). anyway you like.)
    if not m:
        raise click.BadParameter("'" + value + "' is not a range of number. Expected forms like '0-5' or '2'.", 
                           ctx=ctx, param=param)
    start = int(m.group(1), 10)
    end = int(m.group(2) or m.group(1), 10)
    return list(range(start, end + 1))



def str2bool(ctx, param, value):
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if not isinstance(value, str):
        raise click.BadParameter("Boolean value expected.")

    v = value.lower()
    if v in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v in ('no', 'false', 'f', 'n', '0'):
        return False
    raise click.BadParameter("Boolean value expected.")

def multi_str2bool(ctx, param, value):
    return tuple(str2bool(ctx, param, v) for v in value)

def str2none(ctx, param, value):
    if value is None:
        return None
    v = value.lower()
    if v in ('none', 'null'):
        return None
    return value

def multi_str2none(ctx, param, value):
    return tuple(str2none(ctx, param, v) for v in value)


class OrderedGroup(click.Group):
    def list_commands(self, ctx):
        return list(self.commands)
    

def _add_obs(ctx, param, value, obs_name):
    if not value:
        return
    current = list(ctx.params.get("obs") or ())
    current.append(obs_name)
    ctx.params["obs"] = tuple(current)

