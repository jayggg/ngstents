# Spacetime Tents in NGSolve

This `ngstents` package is an c++ extension of the
[NGSolve](https://ngsolve.org) finite element library designed to ease
experimentation with solvers based on spacetime tents for hyperbolic
systems. The python front-end allows new equations (linear or
nonlinear conservation laws) and schemes to be easily defined using
their fluxes and numerical fluxes in a few lines of code.


## Build

Ensure that NGSolve is installed.  After cloning this repository, compile 
the c++ code here and install:

* `cd src`
* `mkdir build && cd build`
* `cmake -DNGSolve_DIR=<Path2YourNGSolveInstall>  ..`
* `make install`

## Check

You can check if your installation is working correctly using the test
suite provided in the 'tests' folder. For example, if you have
`pytest` installed, move to the `tests` folder and use pytest:

* `pytest .`

