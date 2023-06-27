__all__ = ["TentSlab", "Tent", "conslaw", "utils"]

from ._pytents import TentSlab, Tent
from . import conslaw
from . import utils

from .utils._drawtents import DrawPitchedTentsPlt

TentSlab.DrawPitchedTentsPlt = DrawPitchedTentsPlt
del DrawPitchedTentsPlt
