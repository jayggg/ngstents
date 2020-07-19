## NGSTents

This project is based heavily on an earlier 'ConservationLaw' project by Christoph Wintersteiger, Joachim Schöberl and Jay Gopalakrishnan.  The near term goals of this project are to make the code more flexible, for example, to allow new conservation laws, defined by their fluxes and numerical fluxes, to be defined in Python, to allow experimentation with different tent-pitching methods different propagation methods.  Another goal is to implement so-called 'super-tents', each of which will be pitched on a portion of the spatial mesh and solved by an MPI process.


### Installation

NGSolve is required

After cloning this repository, do the following:

* `cd cpp`
* `mkdir build && cd build`
* `cmake ..`
* `make -j4 install`

### Testing

The 'tests' directory contains tests that can be run using pytest.

### Sample Simulations

Currently the only conservation law implemented is the Burgers equation.  The 'burgers' directory contains simple Python scripts for one and two spatial dimensions.

### 1D Visualization

A simple Matplotlib-based code provides visualization of the slab in this case.

### 2D Visualization

Several methods are available for visualizing tent-pitched slabs in two spatial dimensions plus time:

* `tentswebgui` allows visualization of a tent-pitched slab in a browser directly or within a Jupyter notebook.
* tentsngsgui is a plugin to NGSolve/ngsgui which adds a tent scene for visualizing a tent-pitched slab
* The tent slab geometry can also be exported to a VTK file for viewing with Paraview.

### 3D Visualization

The tentsngsgui currently supports a tent-pitched slab on a 3D spatial mesh, but gives very incomplete visualization of the tents and layers.

