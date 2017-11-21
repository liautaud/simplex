"""
A collection of pivots for the Simplex algorithm.
"""
__author__ = 'Romain Liautaud'
__email__ = 'romain.liautaud@ens-lyon.fr'


from random import choice


def pivot_interactive(tb):
    entering_candidates = ', '.join(
        ['x_' + str(i + 1) for i in tb.get_entering_candidates()])
    print('Entering variable candidates:', entering_candidates)
    entering = int(input('Choose the entering variable: ')) - 1

    leaving_candidates = ', '.join(
        ['x_' + str(i + 1) for i in tb.get_leaving_candidates(entering)])
    print('Leaving variable candidates:', leaving_candidates)
    leaving = int(input('Choose the leaving variable: ')) - 1
    print()

    return entering, leaving


def pivot_maximum(tb):
    entering_pairs = [(tb.c[0, j], j) for j in tb.get_entering_candidates()]
    entering_pairs.sort(reverse=True)

    entering = entering_pairs[0][1]
    leaving = tb.get_leaving_candidates(entering)[0]

    return entering, leaving


def pivot_bland(tb):
    entering = tb.get_entering_candidates()[0]
    leaving = tb.get_leaving_candidates(entering)[0]

    return entering, leaving


def pivot_random(tb):
    entering = choice(tb.get_entering_candidates())
    leaving = choice(tb.get_leaving_candidates(entering))

    return entering, leaving


rules = {
    'interactive': pivot_interactive,
    'maximum': pivot_maximum,
    'bland': pivot_bland,
    'random': pivot_random}
