# Spacetime Tents in NGSolve

This `ngstents` package is a c++ extension of the
[NGSolve](https://ngsolve.org) finite element library, designed to ease
experimentation with solvers based on spacetime tents for hyperbolic
systems. A python front-end allows new equations (linear or
nonlinear conservation laws) to be solved by easily defining 
required fluxes and numerical fluxes in a few lines of code.


## Build

Ensure that NGSolve is built.  After cloning this repository, compile 
the c++ code here and install:

* `cd src`
* `mkdir build && cd build`
* `cmake -DNGSolve_DIR=<Path2YourNGSolveInstallCMake>  ../src`
* `make install`

## Use

The `demo` folder contains many examples. 

## Check

New code changes can be cross-checked against a test suite provided in
the 'tests' folder. For example, if you have `pytest` installed, move
to the `tests` folder and use pytest:

* `pytest .`

## Organization

* `src`: c++ source code
* `py`: python packaging files
* `doc`: tutorials and explanations
* `demo`: examples of applications
* `tests`: test suite

