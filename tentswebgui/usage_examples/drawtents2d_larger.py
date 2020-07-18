from netgen.geom2d import unit_square
from ngsolve import Mesh
from tentswebgui import Draw
from ngstents import TentSlab

mesh = Mesh(unit_square.GenerateMesh(maxh=0.05))
print("Fetching tent data")
tps = TentSlab(mesh, dt=0.05, c=5, heapsize=10000000)
print("Generating tents_larger.html")
Draw(tps, 'tents_larger.html')
