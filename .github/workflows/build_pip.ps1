$ErrorActionPreference = "Stop"

if (test-path _skbuild) {
    cmd.exe /c rd /s /q _skbuild
}
if (test-path dist) {
    cmd.exe /c rd /s /q dist
}
if (test-path ..\venv_ngs) {
    cmd.exe /c rd /s /q ..\venv_ngs
}


#$env:NETGEN_CCACHE = 1


#$py=$args[0]

#& $py\python.exe -m venv ..\venv_ngs
#..\venv_ngs\scripts\Activate.ps1
#$env:PATH += ";$env:CI_PROJECT_DIR\venv\bin"
$env:PYDIR = "$env:Python3_ROOT_DIR"
Get-ChildItem $env:PYDIR\include
python --version

pip3 install scikit-build wheel numpy twine mkl-devel==2022.* mkl==2022.*
pip3 install ngsolve


$env:NGSolve_DIR = "$env:Python3_ROOT_DIR\lib\site-packages\ngsolve\cmake"
$env:Netgen_DIR = "$env:Python3_ROOT_DIR\lib\site-packages\netgen\cmake"
$env:LIB = "$env:LIB;$env:Python3_ROOT_DIR\libs"
$env:LIBPATH = "$env:LIBPATH;$env:Python3_ROOT_DIR\libs"
Get-ChildItem $env:Python3_ROOT_DIR\libs
#Get-ChildItem $env:NGSolve_DIR 

Set-Location ../..
Get-ChildItem
python setup.py bdist_wheel -G"Visual Studio 16 2019" -d wheelhouse

#python -m twine upload dist\*.whl
