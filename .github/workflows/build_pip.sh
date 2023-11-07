#! /bin/bash
cd /workspace/ngstents
set -e
yum -y update
yum -y install ninja-build fontconfig-devel tk-devel tcl-devel libXmu-devel mesa-libGLU-devel ccache tree
tree 

py=/opt/python/cp39-cp39/bin/python
cd .github/workflows/ && $py fix_auditwheel_policy.py && cd ../..

rm -rf wheelhouse
mkdir wheelhouse

git config --global --add safe.directory '*'

export ORIGINAL_PATH="$PATH"

for pyversion in 38 39 310 311
do
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    export PATH="$ORIGINAL_PATH:$PYDIR"

    rm -rf _skbuild
    $PYDIR/pip install pytest-check numpy wheel scikit-build mkl==2023.* mkl-devel==2023.* setuptools
    $PYDIR/pip install ngsolve 

    $PYDIR/pip wheel -vvv .
    rename linux_ manylinux_2_17_x86_64.manylinux2014_ ngstents*.whl
    mv ngstents*.whl wheelhouse/
    rm -rf *.whl
    $PYDIR/pip uninstall -y ngsolve netgen-mesher setuptools pytest-check numpy wheel scikit-build mkl mkl-devel 
done
