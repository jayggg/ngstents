from ngsolve import Mesh, Draw
from netgen.csg import unit_cube
from ngstents import TentSlab

mesh = Mesh(unit_cube.GenerateMesh(maxh=0.4))
tps = TentSlab(mesh, 0.2, 8)

Draw(mesh)
Draw(tps)

