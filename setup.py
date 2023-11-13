#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import os
import sys
# import netgen.version
import site

from skbuild import setup
import skbuild.cmaker
from subprocess import check_output
from distutils.sysconfig import get_python_lib
import pkg_resources

from os.path import dirname, isdir, join
import os
import re
import subprocess


def is_canonical(version):
    return re.match(r'^([1-9][0-9]*!)?(0|[1-9][0-9]*)(\.(0|[1-9][0-9]*))*((a|b|rc)(0|[1-9][0-9]*))?(\.post(0|[1-9][0-9]*))?(\.dev(0|[1-9][0-9]*))?$', version) is not None

def get_version():
    """
    Gets the current version number.
    If in a git repository, it is the current git tag.
    Otherwise it is the one contained in the PKG-INFO file.
    """
    version_re = re.compile('^Version: (.+)$', re.M)
    d = dirname(__file__)

    if isdir(join(d, '.git')):
        cmd = 'git describe --always --tags --match v[0-9]*'.split()
        try:
            version = subprocess.check_output(cmd).decode().strip()
        except subprocess.CalledProcessError:
            print('Unable to get version number from git tags')
            exit(1)
        if '-' in version:
            version = '.post'.join(version.split('-')[:2])
        with open(os.devnull, 'w') as fd_devnull:
            subprocess.call(['git', 'status'],
                            stdout=fd_devnull, stderr=fd_devnull)
        if not is_canonical(version[1:]):
            print('Unable to get version number from git tags, using dummy version v0.0.0')
            version="v0.0.0"
        version = version[1:]
    else:
        with open(join(d, 'PKG-INFO')) as f:
            version = version_re.search(f.read()).group(1)
    return version

ngsolve_version = '6.2.2301'
install_requires = [ 'ngsolve >= '+ngsolve_version ] #, 'ngstents >= 0.0.2' ]

py_install_dir = get_python_lib(1,0,'').replace('\\','/')
print("python install dir:",py_install_dir)

_cmake_args = ['-DCMAKE_CXX_COMPILER=ngscxx']

packages=["ngstents"]

if 'darwin' in sys.platform:
    _cmake_args += ['-DPY_INSTALL_DIR='+py_install_dir]
    pass
elif 'linux' in sys.platform:
    install_requires.append('mkl')
elif 'win' in sys.platform:
    install_requires.append('mkl')

cmake_prefix_path = ""
if 'PYDIR' in os.environ and 'win' not in sys.platform:
    print("PYDIR",os.environ["PYDIR"])
    cmake_prefix_path += os.environ["PYDIR"]
    _cmake_args += [f'-DPYTHON_EXECUTABLE={os.environ["PYDIR"]}/python3']
    _cmake_args += [f'-DPYTHON_LIBRARY={os.environ["PYDIR"]}/../lib']
    _cmake_args += [f'-DPYTHON_INCLUDE_DIR={os.environ["PYDIR"]}/../include']
if 'NGSolve_DIR' in os.environ:
    _cmake_args += [f'-DNGSolve_DIR={os.environ["NGSolve_DIR"]}']
if 'Netgen_DIR' in os.environ:
    _cmake_args += [f'-DNetgen_DIR={os.environ["Netgen_DIR"]}']
if 'CMAKE_PREFIX_PATH' in os.environ:
    cmake_prefix_path += os.environ["CMAKE_PREFIX_PATH"]

_cmake_args += ['-DCMAKE_PREFIX_PATH='+cmake_prefix_path]

setup(
    name='ngstents',
    version=str(get_version()),
    author='Jay Gopalakrishnan',
    author_email='gjay@pdx.edu',
    description='Spacetime tent facilities for solving hyperbolic equations',
    long_description='This ngstents package is a c++ extension of the NGSolve finite element library, designed to ease experimentation with solvers based on spacetime tents for hyperbolic systems. A python front-end allows new equations (linear or nonlinear conservation laws) to be solved by easily defining required fluxes and numerical fluxes in a few lines of code.',
    url="https://github.com/jayggg/ngstents",
    license="LGPL2.1",
    install_requires=install_requires,
    packages=packages,
    package_dir={"ngstents": "src"},
    cmake_args=_cmake_args,
    cmake_source_dir='src',
)
