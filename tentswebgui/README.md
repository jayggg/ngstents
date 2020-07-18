##  TentsWebGUI

A GUI for visualizing space-time tents in a web browser, jupyter notebook or jupyterlab context.

### Prerequisites

* Python 3.6 or above
* nodejs and npm (node package manager)
* ipywidgets and widgets notebook extension
`pip3 install --user ipywidgets`
`jupyter nbextension install --py widgetsnbextension`
jupyter nbextension enable --py widgetsnbextension`

### Installation

1. Build the Javascript outputs including 'standalone.js', the notebook widget and the main Python file 'tentswebgui.py' from 'standalone.js' and tentswebgui_template.py.

`./build.sh`

2. Install this package using pip.  For example, in the current directory where 'setup.py' is located, run something like this, depending on your Python setup:

`pip3 install --user .`

Note that the `.` is important; it indicates the current directory.

3. Install the jupyter and jupyterlab extensions

`jupyter nbextension install --user --py tentswebgui`
`jupyter nbextension enable --user --py tentswebgui`


### Testing that everything works

1. In usage_examples, run 'drawtents2d.py' by `python3 drawtents2d.py`.  This should generate a file 'tents.py'.  Open this in a browser, e.g. `firefox drawtents2d.py'

2. You should see a layer of red tetrahedrons viewed from above and some controls in the upper right.  Select to display by layers or tents.  Adjust the sliders by clicking or dragging or typing values in the textboxes.  Use mouse or touch to rotate, zoom and translate the tent scene.

3. To test a larger tent scene on your system, run 'drawtents2d_larger' to generate a file 'tents_larger.py'.   This should open in your browser.  Play with the controls as in the other example.

4. Start the Jupyter notebook server: `jupyter notebook` then open and run Tents2D.ipynb and Burgers2D.ipynb.
