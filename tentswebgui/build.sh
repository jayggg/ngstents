#!/usr/bin/env bash

cd js
# complile Javascript to create standalone.js and .map file
npm run build
cd ..

cp -r package/nbextension src/tentswebgui
cp -r package/labextension src/tentswebgui
cp js/extension.js src/tentswebgui/nbextension/static

# insert compiled Javascript into a copy of webgui_template.py to create
# webgui.py file
echo "creating tentswebgui.py"

python3 build.py . js

cp tentswebgui.py src/tentswebgui
