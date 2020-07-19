##  TentsWebGUI

A GUI for visualizing space-time tents in a web browser, jupyter notebook or (in theory) jupyterlab context.

### Prerequisites

* Python 3.6 or above
* Jupyter
* nodejs and npm (node package manager).  On Mac OS, `brew install node` should work. On Linux there are several options.
* ipywidgets and widgets notebook extension
`pip3 install --user ipywidgets`
`jupyter nbextension install --py widgetsnbextension`
jupyter nbextension enable --py widgetsnbextension`

If you want to try to use Jupyter lab, you need this, I think:

`jupyter labextension install @jupyter-widgets/jupyterlab-manager`

But I haven't been able to get either tentswebgui or NGSolve's webgui working with `jupyter lab` yet.

### Installation

1. Build the Javascript outputs including 'standalone.js', the notebook widget and the main Python file 'tentswebgui.py' from 'standalone.js' and tentswebgui_template.py.

`./build.sh`

You will likely see some warnings, but the build should complete without errors.

2. Install this package using pip.  For example, in the current directory where 'setup.py' is located, run something like this, depending on your Python setup:

`pip3 install --user .`

Note that the `.` is important; it indicates the current directory.

3. Install the jupyter and jupyterlab extensions

`jupyter nbextension install --user --py tentswebgui`
`jupyter nbextension enable --user --py tentswebgui`


### Testing that this project works

1. In the 'demo' folder, run 'drawtents2d.py' by `python3 drawtents2d.py`.  This should generate a file 'tents.html'.  Open this in a browser, e.g. `firefox tents.html` or on Mac OS `open tents.html`.

2. You should see a layer of red tetrahedrons viewed from above and some controls in the upper right.  Select to display by layers or tents.  Adjust the sliders by clicking or dragging or typing values in the textboxes.  Use mouse or touch to rotate, zoom and translate the tent scene.

3. To test a larger tent scene on your system, run 'drawtents2d_larger.py' to generate a file 'tents_larger.html'.   This should open in your browser.  Adjust the controls as in the other example.

4. Start the Jupyter notebook server: `jupyter notebook` then open and run Tents2D.ipynb and Burgers2D.ipynb in the demo folder.  The Burgers2D notebook may take some time to execute.

5. As for the jupyterlab plugin, you can try this in ipython:
``` 
from tentswebgui import *
howtoInstallJupyterLabextension()
```
The latter function generates a command you can try on your system.  It didn't work for me, but it gave a link to an error log.   The error messages seemed to indicate that it neededed npm install to be run again inside the labextension directory, which really seems a bit excessive.
