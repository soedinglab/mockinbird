import argparse
import sys

from mockinbird import __version__
from mockinbird.scripts.preprocess import register_arguments as preprocess_register_args
from mockinbird.scripts.preprocess import run as preprocess_main
from mockinbird.scripts.postprocess import register_arguments as postprocess_register_args
from mockinbird.scripts.postprocess import run as postprocess_main


def create_parser():
    parser = argparse.ArgumentParser(
        'mockinbird',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparser = parser.add_subparsers(title='subcommands')
    for mod_cls in StammpModule.__subclasses__():
        mod_cls(subparser)
    parser.add_argument('--version', action='version',
                        version=__version__)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    if not sys.argv[1:]:
        args = parser.parse_args(['--help'])
    else:
        args = parser.parse_args()

    args.subcommand_func(args)


class StammpModule(object):

    subcommand = None
    aliases = []

    def __init__(self, parser, description=None, help=None):
        subcommand_parser = parser.add_parser(
            self.subcommand,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            description=description,
            help=help,
            aliases=self.aliases
        )

        subcommand_parser.set_defaults(subcommand_func=self)
        self.subcommand_parser = subcommand_parser

    def __call__(self, args):
        pass


class PreprocessModule(StammpModule):
    subcommand = 'preprocess'

    def __init__(self, parser):
        help_msg = 'run preprocessing pipeline'
        description = 'start preprocessing pipeline using a config script'
        super().__init__(parser, help=help_msg, description=description)
        scp = self.subcommand_parser
        preprocess_register_args(scp)

    def __call__(self, args):
        preprocess_main(args)


class PostprocessModule(StammpModule):
    subcommand = 'postprocess'

    def __init__(self, parser):
        help_msg = 'run postprocessing pipeline'
        description = 'start postprocessing pipeline using a config script'
        super().__init__(parser, help=help_msg, description=description)
        scp = self.subcommand_parser
        postprocess_register_args(scp)

    def __call__(self, args):
        postprocess_main(args)


if __name__ == '__main__':
    main()
