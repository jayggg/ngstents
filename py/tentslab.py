from ngsolve import Mesh
from netgen.geom2d import unit_square
from tents import TPS2

mesh = Mesh(unit_square.GenerateMesh(maxh=0.5))

t = TPS2(mesh, 0.1, 1)
print(t.GetNTents())
