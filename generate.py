#!/usr/bin/env python3
"""
A simple generator of linear programming problems.
"""
__author__ = 'Romain Liautaud'
__email__ = 'romain.liautaud@ens-lyon.fr'


import argparse
import random


def generate(n, m, width, density):
    """Print a random linear problem."""

    def generate(k):
        """Generate k random coefficients in a row."""
        coefficients = []

        for i in range(k):
            if random.random() <= density:
                c = 0
            else:
                c = random.randint(-width, width)

            coefficients.append(c)

        return coefficients

    print(n)
    print(m)
    print(' '.join(map(str, generate(n))))
    print(' '.join(map(str, generate(m))))

    for i in range(m):
        print(' '.join(map(str, generate(n))))


def main():
    """Parse the command-line arguments and generate."""

    parser = argparse.ArgumentParser(
        description='Generates linear programming problems.')

    parser.add_argument(
        'n', metavar='n', type=int, help='The number of variables.')

    parser.add_argument(
        'm', metavar='m', type=int, help='The number of constraints.')

    parser.add_argument(
        '--width', type=int, default=100,
        help='The range of the randomly generated coefficients.')

    parser.add_argument(
        '--density', type=float, default=.2,
        help='The probability that a coefficient is equal to 0.')

    args = parser.parse_args()
    generate(args.n, args.m, args.width, args.density)


if __name__ == '__main__':
    main()
