import argparse
import sys
import time


class ArgumentParserError(Exception): pass

class ThrowingArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise ArgumentParserError(message)

parser = ThrowingArgumentParser()

parser.add_argument("a", type=float, help="Premier nombre")
parser.add_argument("b", type=float, help="Deuxieme nombre")
parser.add_argument("tolerance", type=float, help="Tolerance")


try:
	args=parser.parse_args()
except ArgumentParserError, e:
	sys.exit(1)


a=args.a
b=args.b
tolerance=args.tolerance

if (abs(a-b)>=tolerance):
	sys.exit(1)
else:
	sys.exit(0)


sys.exit(1)

