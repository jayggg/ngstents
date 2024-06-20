# Spacetime Tents in NGSolve

This NGS-Tents (`ngstents`) package is a c++ extension of the
[NGSolve](https://ngsolve.org) finite element library, designed to ease
experimentation with solvers based on spacetime tents for hyperbolic
systems. A python front-end allows new equations (linear or
nonlinear conservation laws) to be solved by easily defining 
required fluxes and numerical fluxes in a few lines of code.



## Install



### Build using pip

On a computer with a build system (compiler etc), you can install  NGS-Tents by 

``` sh
python3 -m pip install git+https://github.com/jayggg/ngstents.git
```

If you do not have the dependency [`ngsolve`](https://ngsolve.org)
installed, this command will attempt to install ngsolve first, before
proceeding to install `ngstents`. If you prefer to use your own
existing install of ngsolve, please use the `--no-build-isolation`
argument:

```sh
python3 -m pip install --no-build-isolation git+https://github.com/jayggg/ngstents.git
```

### Build using CMake

If you built ngsolve from source, you can also build and install
`ngstents` in the traditional manner. After cloning this repository,
compile the c++ code here and install as follows:

```sh
cd src
mkdir build && cd build
cmake -DNGSolve_DIR=<Path2YourNGSolveInstallCMake>  ..
make install
```

(Often CMake is able to correctly detect the path to your NGSolve installation in which case you do not have to specify the `NGSolve_DIR` variable.)


### Binary install 

If you do not have a compiler, then you can install NGS-Tents using
a binary installer. To do so on linux, mac, and windows (with python>= 3.9),
use the following command.

```sh
pip install --pre ngstents
```

	



## Use

To start using this code's python interface, import the module after installation:

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

``` sh
pip install -r requirements.txt
sphinx-build -b html . _build_docs
```

You can then open the docs from `ngstents/docs/_build_docs/INDEX.html`. 

Alternately, you can build a jupyter book of the documentation by
navigating to the `ngstents/docs` folder and issuing 

``` sh
jupyter-book build .
```

which creates documentation in `ngstents/docs/_build/INDEX.html`. 


## Check

If you extend the code, your new code changes can be cross-checked
against a test suite provided in the 'tests' folder. For example, if
you have `pytest` installed, move to the `ngstents/tests` folder and use
pytest:

``` sh
pytest .
```


## Organization

* `src`: c++ source code
* `py`: python packaging files
* `doc`: tutorials and explanations
* `demo`: examples of applications
* `tests`: test suite

