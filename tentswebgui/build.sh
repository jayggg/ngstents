#!/usr/bin/env bash

cd js
# compile Javascript to create standalone.js and .map file
# and files for Jupyter
npm run build
cp -r nbextension ../src/tentswebgui
cp -r labextension ../src/tentswebgui
rm -rf nbextension labextension

cp extension.js ../src/tentswebgui/nbextension/static

cd ..
# insert compiled Javascript into a copy of webgui_template.py to create
# webgui.py file
echo "creating tentswebgui.py"

python3 py/build.py py js 

mv py/tentswebgui.py src/tentswebgui
