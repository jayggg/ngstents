from ngsolve import Mesh
from netgen.geom2d import unit_square
from ngstents import TentSlab


# Construct a mesh of tents

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
tps = TentSlab(mesh, dt=0.1, c=10)
print('\nTent pitched spacetime slab made:')
print('   Number of tents:', tps.GetNTents())
print('   Maximal tent slope:', tps.MaxSlope())

# Query tents

n = 100
t = tps.GetTent(n)
print('Details of Tent #%d:' % n)
print('  Pitching (central) vertex number:', t.vertex)
print('  Neighbor vertex numbers:', list(t.nbv))
print('  Tent element numbers:', list(t.els))
print('  Neighbor vertex heights:', list(t.nbtime))
