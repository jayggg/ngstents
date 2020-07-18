from ngsolve import Mesh, Draw
from netgen.csg import unit_cube
from ngstents import TentSlab

mesh = Mesh(unit_cube.GenerateMesh(maxh=0.4))
tps = TentSlab(mesh, 0.2, 8)

Draw(mesh)
Draw(tps)

# Leave this out for now, AdvancingFront method doesn't exist.
# But it might be useful to add this, especially in 3D case.

# Draw(cl.Tau())
# input('start')
# for i in range(cl.GetNTents()):
#     cl.AdvancingFront(i)
# input('step')
# Draw(cl.Tau())
# print(min(cl.Tau().vec),max(cl.Tau().vec))
