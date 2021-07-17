# Spacetime Tents in NGSolve: `ngstents`

`ngstents` is an c++ extension of the [NGSolve](https://ngsolve.org)
finite element library designed to ease experimentation with
tent-based methods for hyperbolic systems. The python front-end allows
new equations (linear or nonlinear conservation laws) and schemes to
be easily defined using their fluxes and numerical fluxes in a few
lines of code.


## Installation

NGSolve is required

After cloning this repository, do the following:

* `cd cpp`
* `mkdir build && cd build`
* `cmake ..`
* `make -j4 install`

Check if your installation is working by running the code in 
the 'tests' fodler (say, by using an automated tester like
pytest).

## Visualization

### 1D + time 

A simple Matplotlib-based code provides visualization of  a (2D) spacetime 
slab in this case.

### 2D + time 

Several methods are available for visualizing tent-pitched slabs in two spatial dimensions plus time:

* `tentswebgui` allows visualization of a tent-pitched slab in a browser directly or within a Jupyter notebook.
* `tentsngsgui` is a plugin to NGSolve/ngsgui which adds a tent scene for visualizing a tent-pitched slab
* The tent slab geometry can also be exported to a VTK file for viewing with Paraview.

### 3D + time

`tentsngsgui` currently supports a tent-pitched slab on a 3D spatial mesh, but gives very incomplete visualization of the tents and layers.

For more details on `tentsngsgui` and `tentswebgui`, see the README.md files in their directories.
