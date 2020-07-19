## TentsNgsGUI

This is a plugin to the [ngsgui](https://github.com/NGSolve/ngsgui) GUI for NGSolve.  It renders using OpenGL.  The documentation for ngsgui is a bit out of date.  To install it you only need to do: `pip3 install ngsgui`  This will install pyside 5.14, which works with Python 3.8 on Mac OS or Linux.

To ensure that `ngsgui` has installed successfully, you can navigate to the NGSolve py_tutorial in the source tree and try `ngsolve poisson.py`

To install this plugin, simply do `pip3 install --user .`

The `.` indicates the current directory (where setup.py is contained).

Once this is installed, if you construct a TentSlab as e.g. 'ts', add a line in your python script `Draw(ts)`, then run your script like `ngsolve myscript.py`, you should see a TentScene with your tent slab along with scenes for any other NGSolve objects you have called `Draw` on.

There are also demo scripts in the demo folder that can be run, via e.g. `ngsolve drawtents2d.py`.


