from netgen.geom2d import unit_square
from ngsolve import Mesh
from tentswebgui import Draw
from ngstents import TentSlab

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
print("Fetching tent data")
tps = TentSlab(mesh, dt=0.1, c=5)
print("Generating tents.html")
Draw(tps, 'tents.html')
