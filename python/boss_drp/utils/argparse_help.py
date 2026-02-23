
import argparse
import sys

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
