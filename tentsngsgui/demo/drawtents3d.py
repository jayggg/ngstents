from ngsolve import Mesh, Draw
from netgen.csg import unit_cube
from ngstents import TentSlab

mesh = Mesh(unit_cube.GenerateMesh(maxh=0.4))

tentslab = TentSlab(mesh, "edge")
c = 1
tentslab.SetMaxWavespeed(c)
tentslab.PitchTents(dt=0.1, local_ct=True, global_ct=1)

Draw(mesh)
Draw(tentslab)

