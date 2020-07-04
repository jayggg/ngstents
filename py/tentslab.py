from ngsolve import Mesh
from netgen.geom2d import unit_square
from tents import TPS2

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))

t = TPS2(mesh, dt=0.1, c=10)
print('\nTent pitched spacetime slab made:')
print('   Number of tents:', t.GetNTents())
print('   Maximal tent slope:', t.MaxSlope())
