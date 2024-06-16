__all__ = ["TentSlab", "Tent", "conslaw", "utils"]

import ngsolve
from ._pytents import TentSlab, Tent
from . import conslaw
from . import utils

from .utils._drawtents import DrawPitchedTentsPlt
from .utils._drawtents2d import DrawPitchedTents

TentSlab.DrawPitchedTentsPlt = DrawPitchedTentsPlt
TentSlab.DrawPitchedTents = DrawPitchedTents
del DrawPitchedTentsPlt, DrawPitchedTents
