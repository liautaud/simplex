#!/usr/bin/env python3
"""
A Python implementation of the Simplex Algorithm for the OA course.

It supports the following pivot rules:
- Maximum Coefficient Rule (using --rule=pivot)
- Bland's Rule (using --rule=bland)
- Random pivot (using --rule=random)
"""
__author__ = 'Romain Liautaud'
__email__ = 'romain.liautaud@ens-lyon.fr'
__license__ = 'GPL'


import argparse
import sympy

from sympy.matrices import Matrix


class Program:
    def __init__(self, n, m, c, b, a):
        """
        Initialize a linear program.

        n: Number of variables of the program.
        m: Number of constraints of the program.
        c: List of coefficients for the objective function.
        b: List of coefficients for the right-hand side.
        a: List of list of coefficients for the left-hand side.
        """
        self.n = n
        self.m = m
        self.c = Matrix(c)
        self.b = Matrix(b)
        self.a = Matrix(a)

    def print(self):
        """
        Pretty-print the linear program.
        """
        var_padding = len(str(self.n))
        coe_padding = max(
            *[len(str(k)) for k in self.c],
            *[len(str(k)) for k in self.b],
            *[len(str(k)) for k in self.a])
        padding = 2 + var_padding + coe_padding

        def align(row):
            s = ''
            for i, k in enumerate(row):
                if k == 0:
                    s += ' ' * (padding + 1)
                else:
                    t = str(k) + '_x' + str(i + 1)
                    if not t.startswith('-'):
                        t = '+' + t
                    s += t.ljust(padding + 1)
            return s

        print('Maximize  ' + align(self.c))

        for i in range(self.m):
            if i == 0:
                row_start = 'Such that '
            else:
                row_start = ' ' * 10

            print(row_start +
                  align(self.a.row(i)) +
                  '<= ' + str(self.b[i]))

        var_names = ', '.join(['x_' + str(i + 1) for i in range(self.n)])
        print((' ' * 10) + var_names + ' are non-negative')

    @staticmethod
    def parse(data):
        """
        Parse a linear program from an input stream.

        data: Input stream, as returned by open().
        """
        n = int(data.readline())
        m = int(data.readline())

        lines = []
        for i in range(m + 2):
            line = data.readline().strip().split(' ')
            line = list(map(sympy.Rational, line))
            lines.append(line)

        if len(lines[0]) != n:
            raise Exception('Not enough coefficients for the obj. function.')
        if len(lines[1]) != m:
            raise Exception('Not enough coefficients for the RHS.')
        if any([len(lines[i]) != n for i in range(2, m + 2)]):
            raise Exception('Not enough coefficients for the LHS.')

        return Program(n, m, lines[0], lines[1], lines[2:])


class Tableau:
    pass


def main():
    """
    Parse the command-line arguments and run the algorithm.
    """
    parser = argparse.ArgumentParser(
        description='Solves a linear problem using the Simplex Algorithm.')

    parser.add_argument(
        'input', metavar='I', type=open,
        help='The input file in the right format.')

    parser.add_argument(
        '--rule', metavar='R', choices=['pivot', 'bland', 'random'],
        help='The pivot rule to use in the algorithm.', default='pivot')

    parser.set_defaults(verbose=False)
    parser.add_argument(
        '--verbose', '-v', action='store_true',
        help='Whether to enable verbose output.')

    args = parser.parse_args()
    # print(args)
    Program.parse(args.input).print()


if __name__ == '__main__':
    main()
