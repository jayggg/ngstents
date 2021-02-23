from ._pytents import TentSlab, Tent
from . import conslaw
from . import utils

try:
    from ._drawtents import DrawPitchedTentsPlt, Draw3DTentPlt
    TentSlab.DrawPitchedTentsPlt = DrawPitchedTentsPlt
    TentSlab.Draw3DTentPlt = Draw3DTentPlt
    del DrawPitchedTentsPlt
    del Draw3DTentPlt
except ModuleNotFoundError:
    pass
