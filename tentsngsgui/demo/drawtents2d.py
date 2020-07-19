from ngsolve import Mesh, Draw
from netgen.geom2d import unit_square
from ngstents import TentSlab

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

tps = TentSlab(mesh, 0.2, 1/0.18)

Draw(mesh)
# Draw(tps)

# Set this aside for now.  AdvancingFront method doesn't exist.
# Draw(cl.Tau())
# input('start')
# for i in range(cl.GetNTents()):
#     cl.AdvancingFront(i)
    # input('step')
# Draw(cl.Tau())
# print(min(cl.Tau().vec),max(cl.Tau().vec))
