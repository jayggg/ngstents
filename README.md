# Spacetime Tents in NGSolve

This `ngstents` package is a c++ extension of the
[NGSolve](https://ngsolve.org) finite element library, designed to ease
experimentation with solvers based on spacetime tents for hyperbolic
systems. A python front-end allows new equations (linear or
nonlinear conservation laws) to be solved by easily defining 
required fluxes and numerical fluxes in a few lines of code.



## Install

Binary installers are available for linux, mac, and windows (with
python >= 3.9).

* `pip install ngstents`


## Build

If you built NGSolve from source, you can build and install `ngstents`.
After cloning this repository, compile  the c++ code here and install:

* `cd src`
* `mkdir build && cd build`
* `cmake -DNGSolve_DIR=<Path2YourNGSolveInstallCMake>  ..`
* `make install`

## Use

To start using this code, import the module after installation:

``` python
import ngstents
```

It contains classes named `Tent` and `TentSlab` used to partition
spacetime into tents. It also has a submodule named `conslaw` which
contains specific pre-programmed hyperbolic equations and also
facilities to define and solve new hyperbolic systems using the python
interface. Use python's help system to query for documentation on each
object.  The [`demo` folder](./demo) contains example scripts. The
[`docs` folder](./docs/INDEX.ipynb) contains more documentation,
explained next.




## Documentation

[Read the docs](https://jayggg.github.io/ngstents/) online. This
documentation is generated from hands-on style tutorial notebooks in
the [`docs` folder](./docs/INDEX.ipynb).

Offline, to build and test the docs locally, navigate to the `ngstents/docs` 
folder and run these:

* `pip install -r requirements.txt`
* `sphinx-build -b html . _build_docs`

You can then open the docs from `ngstents/docs/_build_docs/INDEX.html`. 


## Check

If you extend the code, your new code changes can be cross-checked
against a test suite provided in the 'tests' folder. For example, if
you have `pytest` installed, move to the `ngstents/tests` folder and use
pytest:

* `pytest .`

## Organization

* `src`: c++ source code
* `py`: python packaging files
* `doc`: tutorials and explanations
* `demo`: examples of applications
* `tests`: test suite

