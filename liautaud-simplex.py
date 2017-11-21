#!/usr/bin/env python3
"""
An implementation of the Simplex algorithm for the OA course.

It supports the following pivot rules:
- Interactive pivot (using --rule=interactive)
- Maximum Coefficient Rule (using --rule=maximum)
- Bland's Rule (using --rule=bland)
- Random pivot (using --rule=random)
"""
__author__ = 'Romain Liautaud'
__email__ = 'romain.liautaud@ens-lyon.fr'


import argparse
from enum import Enum
from sympy import Matrix, Rational, nan, oo

from pivots import rules


class State(Enum):
    UNFEASIBLE = 0
    BOUNDED = 1
    UNBOUNDED = 2


class Program:

    def __init__(self, n, m, c, b, a):
        """Initialize a linear program.

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
        """Pretty-print the linear program."""

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
        """Parse a linear program from an input stream.

        data: Input stream, as returned by open().
        """
        n = int(data.readline())
        m = int(data.readline())

        lines = []
        for i in range(m + 2):
            line = data.readline().strip().split(' ')
            line = list(map(Rational, line))
            lines.append(line)

        if len(lines[0]) != n:
            raise Exception('Not enough coefficients for the obj. function.')
        if len(lines[1]) != m:
            raise Exception('Not enough coefficients for the RHS.')
        if any([len(lines[i]) != n for i in range(2, m + 2)]):
            raise Exception('Not enough coefficients for the LHS.')

        return Program(n, m, lines[0], lines[1], lines[2:])


class Tableau:

    def __init__(self, c, basic, nonbasic, artificial):
        """Initialize a tableau.

        c: List of coefficients of the tableau.
        basic: List of basic variables.
        nonbasic: List of non-basic variables.
        artificial: List of artificial variables.
        """
        self.c = Matrix(c)
        self.basic = basic
        self.nonbasic = nonbasic
        self.artificial = artificial

    def get_objective_value(self):
        """Return the value of the objective function."""
        return -1 * self.c[0, -1]

    def set_objective(self, v):
        """Redefine the objective function."""
        self.c[0, :-1] = [v + [0] * (self.c.cols - 1 - len(v))]

        # We start with an objective value of 0, which is what we would get
        # if all variables were non-basic. We will update this value as we
        # substitute basic variables for non-basic ones later.
        self.c[0, -1] = 0

        self.do_pivot_objective()

    def get_entering_candidates(self):
        """Return the variables which can be chosen to enter the basis."""
        return [j for j in self.nonbasic if self.c[0, j] > 0]

    def get_leaving_candidates(self, entering):
        """Return the variables which can be chosen to leave the basis."""
        candidates = []
        minimum_ratio = oo

        # Perform the ratio test on every row, and return the ties.
        for row, i in enumerate(self.basic):
            if self.c[row + 1, entering] > 0:
                ratio = self.c[row + 1, -1] / self.c[row + 1, entering]
                if ratio < minimum_ratio:
                    minimum_ratio = ratio
                    candidates = [i]
                elif ratio == minimum_ratio:
                    candidates.append(i)

        print('Leaving candidates:', candidates)

        return candidates

    def is_optimal(self):
        """Return whether the tableau meets the optimality criterion."""
        return len(self.get_entering_candidates()) == 0

    def is_unbounded(self):
        """Return whether the tableau meets the unbounded criterion."""
        return any([self.c[0, j] > 0 and
                    all([self.c[i, j] <= 0 for i in range(1, self.c.rows)])
                    for j in range(self.c.cols)])

    def remove_artificial(self):
        """Remove all the artificial variables from the tableau.

        Some of those variables might still be in the basis, so we might have
        have to pivot a couple of times to remove them from the basis before
        deleting the corresponding columns.
        """
        available = [i for i in self.nonbasic if i not in self.artificial]

        for j in self.artificial:
            if j in self.basic:
                row = self.basic.index(j)

                for i in available:
                    if self.c[row, i] != 0:
                        self.do_pivot(i, j)
                        available.remove(i)

                self.basic.remove(j)
            else:
                self.nonbasic.remove(j)

            self.c.col_del(j)

    def do_pivot(self, entering, leaving):
        """Apply the pivot transformation on the tableau.

        entering: Variable which will enter the basis.
        leaving: Variable which will leave the basis.
        """
        entering_row = self.nonbasic.index(entering)
        leaving_row = self.basic.index(leaving)

        self.nonbasic[entering_row] = leaving
        self.basic[leaving_row] = entering

        self.gauss_jordan(leaving_row + 1, entering)
        self.do_pivot_objective()

    def gauss_jordan(self, i, j):
        """Apply Gauss-Jordan elimination to the tableau.

        This will transform the j-th column into a column filled with
        zeroes except on the i-th row, which will have coefficient 1.

        Note that this doesn't affect the objective vector: to perform a
        full pivot, you must also call do_pivot_objective() to express
        the objective vector in terms of the non-basic variables only.
        """

        # Start by putting coefficient 1 in c[i, j].
        if self.c[i, j] != 0:
            self.c[i, :] *= Rational(1, self.c[i, j])
        else:
            p = next(a for a in range(self.c.rows)
                     if self.c[a, j] != 0 and a != i)
            self.c[i, :] += self.c[p, :] * Rational(1, self.c[p, j])

        # Then put coefficient 0 in all the other c[b, j].
        for b in range(self.c.rows):
            if b != i:
                self.c[b, :] -= self.c[i, :] * self.c[b, j]

    def do_pivot_objective(self):
        """Express the objective vector in terms of the non-basic variables.

        This is used in the set_objective and do_pivot methods to maintain
        the invariant needed by the Simplex.
        """

        # To do so, we will substitute the basic variables with a linear
        # combination of the non-basic ones and a scalar.
        for row, i in enumerate(self.basic):
            if self.c[0, i] != 0:
                a = self.c[0, i]

                for j in range(self.c.cols - 1):
                    self.c[0, j] -= a * self.c[row + 1, j]

                self.c[0, -1] -= a * self.c[row + 1, -1]

    def print(self):
        """Pretty-print the tableau."""

        padding = max([len(str(k)) for k in self.c])

        def align(row):
            return ' '.join([str(k).ljust(padding) for k in row[:-1]]) +\
                   ' | ' + str(row[-1])

        print(align(self.c.row(0)))
        print('-' * (1 + self.c.cols * (padding + 1)))

        for i in range(1, self.c.rows):
            print(align(self.c.row(i)))

    @staticmethod
    def convert(prog):
        """Create the initial tableau for a given program.

        This tableau is meant for phase I of the simplex, as it contains
        artificial variables which should not appear in the final solution.

        Warning: this method doesn't fill the objective vector. You must do
        it manually using tb.set_objective().
        """
        n = prog.n
        basic = []
        artificial = []

        # Stores the slack variable that was created for each constraint.
        slack = {}

        # Create the basic tableau from the program.
        c = prog.a[:, :]
        c = c.row_insert(0, Matrix([nan] * n).T)

        def make_col(i):
            """Create a vector with zeros everywhere but on the i-th row."""
            return Matrix([[nan] + [0] * (i - 1) + [1] + [0] * (prog.m - i)])

        # Add the slack variables.
        for i in range(prog.m):
            slack[i] = n
            c = c.col_insert(n + 1, make_col(i + 1).T)
            n += 1

        # Make all the right-hand side coefficients positive, and add
        # artificial variables if needed. At the end, we get as many
        # basic variables as there are constraints: they are either
        # the slack variables we created earlier, or the artificial
        # variables we just created.
        for i in range(prog.m):
            if prog.b[i] < 0:
                c[i + 1, :] *= -1
                c = c.col_insert(n + 1, make_col(i + 1).T)
                artificial.append(n)
                basic.append(n)
                n += 1
            else:
                basic.append(slack[i])

        # Add the right-hand side.
        c = c.col_insert(n, Matrix([nan] + list(map(abs, prog.b[:]))))

        # Deduce the non-basic variables.
        nonbasic = [i for i in range(n) if i not in basic]

        return (len(artificial), Tableau(c, basic, nonbasic, artificial))


class Solver:

    def __init__(self, rule, verbose):
        """Initialize a Simplex solver with options.

        rule: The pivot rule to use.
        verbose: Whether the output should be verbose.
        """
        self.rule = rule
        self.verbose = verbose

    def solve(self, prog):
        """Solve the given linear program."""
        a_count, tb = Tableau.convert(prog)

        if a_count > 0:
            # Phase I of the algorithm.
            na_count = tb.c.cols - 1 - a_count
            initial_obj = [0] * na_count + [-1] * a_count
            tb.set_objective(initial_obj)

            if self.verbose:
                print('[Starting phase I]')
                print('The phase I tableau is:')
                print()
                tb.print()
                print()

            state, tb, _ = self.one_phase(tb, False)

            if state != State.BOUNDED:
                print('This linear program is UNFEASIBLE.')
                return

            if tb.get_objective_value() != 0:
                print('This linear program is UNFEASIBLE.')
                return

            tb.remove_artificial()

        # Phase II of the algorithm.
        tb.set_objective(prog.c[:])

        if self.verbose:
            print('[Starting phase II]')
            print('The phase II tableau is:')
            print()
            tb.print()
            print()

        state, tb, pivots = self.one_phase(tb)

        if state == State.BOUNDED:
            solution = self.get_solution(prog, tb)
            solution_str = ', '.join([
                'x_' + str(i + 1) +
                ' = ' + str(solution[i]) for i in solution])

            print('This linear program has a BOUNDED solution.')
            print('One optimal solution is: %s' % solution_str)
            print('The value of the objective for this solution is: %d' %
                  tb.get_objective_value())
            print('The number of pivots is: %d' % pivots)
            print('The pivot rule used: %s' % self.rule)
            return

        if state == State.UNBOUNDED:
            print('This linear program has an UNBOUNDED solution.')
            return

    def one_phase(self, tb, silent=False):
        """Perform one phase of the algorithm, starting on a given tableau.

        silent: Whether to stay quiet no matter the value of verbose.
        """
        pivots = 0

        if not silent and self.verbose:
            print('The initial tableau is:')
            print()
            tb.print()
            print()

            basis = map(lambda i: 'x_' + str(i + 1), tb.basic)
            print('The initial basis is: ' + ', '.join(basis))
            print()

        while True:
            if tb.is_optimal():
                return (State.BOUNDED, tb, pivots)

            if tb.is_unbounded():
                return (State.UNBOUNDED, tb, pivots)

            entering, leaving = rules[self.rule](tb)
            tb.do_pivot(entering, leaving)
            pivots += 1

            if not silent and self.verbose:
                print('The entering variable is: x_%d' % (entering + 1))
                print('The leaving variable is: x_%d' % (leaving + 1))
                print()
                tb.print()
                print()

    def get_solution(self, prog, final_tb):
        """Return the solution for a given program and final tableau."""
        solution = {}

        for i in range(prog.n):
            if i in final_tb.basic:
                index = final_tb.basic.index(i)
                solution[i] = final_tb.c[index + 1, -1]
            else:
                solution[i] = 0

        return solution


def main():
    """Parse the command-line arguments and run the algorithm."""

    parser = argparse.ArgumentParser(
        description='Solves a linear problem using the Simplex algorithm.')

    parser.add_argument(
        'input', metavar='I', type=open,
        help='The input file in the right format.')

    parser.add_argument(
        '--rule', metavar='R', choices=rules.keys(),
        help='The pivot rule to use in the algorithm.', default='bland')

    parser.set_defaults(verbose=False)
    parser.add_argument(
        '--verbose', '-v', action='store_true',
        help='Whether to enable verbose output.')

    args = parser.parse_args()
    program = Program.parse(args.input)
    solver = Solver(args.rule, args.verbose)

    print('The input linear program is:')
    print()
    program.print()
    print()
    solver.solve(program)


if __name__ == '__main__':
    main()
