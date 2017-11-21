# Linear solver - OA course - ENSL 2017

A Python implementation of the Simplex algorithm for the OA course.


## Usage.

```
╒═══════════════════════════════════════╕
| Linear solver - OA course - ENSL 2017 |
╘═══════════════════════════════════════╛

usage: liautaud-simplex.py [-h] [--rule R] [--verbose] I

Solves a linear problem using the Simplex algorithm.

positional arguments:
  I              The input file in the right format.

optional arguments:
  -h, --help     show this help message and exit
  --rule R       The pivot rule to use in the algorithm.
  --verbose, -v  Whether to enable verbose output.
```

When no input file is specified, the input is read from stdin, which allows for:
```
./generate.py {n} {m} | ./liautaud-simplex.py
```


## Installation guide.

Before using the solver for the first time, run `python3 setup.py install --user`
to install the project's dependencies, or install `sympy` yourself using `pip`.


## About this solver.

This solver depends on `argparse` for the parsing of command-line arguments, and
on `sympy` to handle matrices and rational numbers.

It has been tested on problems up to 100 variables, on other problems up to 100
constraints, and on problems with 50 variables and 50 constraint -- even though
this last one took a long while to complete.

See `liautaud-discussion.txt` for remarks and discussion on interesting examples.


## Generating random problems.

To generate random problems, use `./generate.py`, which has the following specs:

```
usage: generate.py [-h] [--width WIDTH] [--density DENSITY] n m

Generates linear programming problems.

positional arguments:
  n                  The number of variables.
  m                  The number of constraints.

optional arguments:
  -h, --help         show this help message and exit
  --width WIDTH      The range of the randomly generated coefficients.
  --density DENSITY  The probability that a coefficient is equal to 0.
```