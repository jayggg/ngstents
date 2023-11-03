$ErrorActionPreference = "Stop"

if (test-path _skbuild) {
    cmd.exe /c rd /s /q _skbuild
}
if (test-path dist) {
    cmd.exe /c rd /s /q dist
}
if (test-path venv_ngs) {
    cmd.exe /c rd /s /q venv_ngs
}

pip3 install scikit-build wheel numpy twine mkl-devel==2022.* mkl==2022.*
pip3 install ngsolve

$env:NGSolve_DIR = "$env:Python3_ROOT_DIR\lib\site-packages\ngsolve\cmake"
$env:Netgen_DIR = "$env:Python3_ROOT_DIR\lib\site-packages\netgen\cmake"

Set-Location ../..
python setup.py bdist_wheel -G"Visual Studio 16 2019" -d wheelhouse
