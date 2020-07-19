## NGSTents

This project is based heavily on an earlier 'ConservationLaw' project by Christoph Wintersteiger, Joachim Schöberl and Jay Gopalakrishnan.  The near term goals of this project are to make the code more flexible, for example, to allow new conservation laws, defined by their fluxes and numerical fluxes, to be defined in Python, to allow experimentation with different tent-pitching methods different propagation methods.

### Installation

NGSolve is required

After cloning this repository, do the following:

* `cd cpp`
* mkdir build && cd build
* cmake ..
* make -j4 install

### Testing

The 'tests' directory contains tests that can be run using pytest.


### Simulations

Currently the only conservation law implmented is the Burgers equation.  The 'burgers' directory contains scripts for one and two spatial dimensions.

