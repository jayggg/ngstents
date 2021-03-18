from ngsolve import Mesh, Draw
from netgen.geom2d import unit_square
from ngstents import TentSlab

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

tentslab = TentSlab(mesh, "edge")
c = 1
tentslab.SetMaxWavespeed(c)
tentslab.PitchTents(dt=0.1, local_ct=True, global_ct=1)

Draw(mesh)
Draw(tentslab)

