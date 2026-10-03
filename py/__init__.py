__all__ = ["TentSlab", "Tent", "conslaw", "utils"]

import os
import sys

import ngsolve

# Keep the core DLL discoverable while importing the extensions on Windows.
if sys.platform.startswith("win"):
    _dll_directory = os.add_dll_directory(os.path.dirname(__file__))
del os, sys

from ._pytents import TentSlab, Tent
from . import conslaw
from . import utils

from .utils._drawtents import DrawPitchedTentsPlt
from .utils._drawtents2d import DrawPitchedTents

TentSlab.DrawPitchedTentsPlt = DrawPitchedTentsPlt
TentSlab.DrawPitchedTents = DrawPitchedTents
del DrawPitchedTentsPlt, DrawPitchedTents
