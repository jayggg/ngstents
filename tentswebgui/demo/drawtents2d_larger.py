from netgen.geom2d import unit_square
from ngsolve import Mesh
from tentswebgui import Draw
from ngstents import TentSlab

mesh = Mesh(unit_square.GenerateMesh(maxh=0.05))
print("Fetching tent data")

# using causality constant
local_ctau = True
global_ctau = 1
wavespeed = 5
dt = 0.05   
ts = TentSlab(mesh, method="edge")
ts.SetMaxWavespeed(wavespeed)
ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())
print("Generating tents.html")
Draw(ts, 'tents.html')
