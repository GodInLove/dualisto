import argparse

import sys

import index
import quant

from version import version

APP_VERSION = (
        '''
    ----------------------------------------------------------------------

    dualisto-v(%s) - an independent demo of dual rna seq data analysis

    ----------------------------------------------------------------------
    ''' % version
)


def prepare(argv):
    parser = argparse.ArgumentParser(description=APP_VERSION,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub_parser = parser.add_subparsers(title="dualisto options")

    # index
    dualisto_index = sub_parser.add_parser("index", help="build a dualisto index")
    dualisto_index.set_defaults(func=index.start)
    index.args_handle(dualisto_index)

    # quant
    dualisto_quant = sub_parser.add_parser("quant", help="run with dualisto")
    dualisto_quant.set_defaults(func=quant.start)
    quant.args_handle(dualisto_quant)

    if len(argv) == 1:
        print(parser.print_help())
        sys.exit(0)
    else:
        all_args = parser.parse_args(argv[1:])
        all_args.func(all_args)


def main():
    prepare(sys.argv)


if __name__ == '__main__':
    main()
