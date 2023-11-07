#! /bin/bash
set -e
cd ../..
rm -rf _skbuild dist venv_ngs

export PYDIR=$Python3_ROOT_DIR/bin

$PYDIR/python3 --version

export PATH=/Applications/CMake.app/Contents/bin:$PATH
export NETGEN_Dir=$PYDIR/../lib/python$1/site-packages/netgen/cmake
export NGSolve_Dir=$PYDIR/../lib/python$1/site-packages/ngsolve/cmake
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$NGSolve_Dir:$NETGEN_Dir
$PYDIR/pip3 install scikit-build wheel

export CMAKE_OSX_ARCHITECTURES='arm64;x86_64'
$PYDIR/pip3 install ngsolve

$PYDIR/python3 setup.py bdist_wheel --plat-name macosx-10.15-universal2 -d wheelhouse
