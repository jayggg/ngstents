from netgen.geom2d import SplineGeometry
from ngsolve import (Mesh, Draw, Redraw, CoefficientFunction, sqrt, x, y,
                     L2, GridFunction, TaskManager)
from ngsolve.internal import visoptions
from ngstents import TentSlab
from ngstents.conslaw import Euler

fine_mesh = False
geom = SplineGeometry()
if(fine_mesh):
    maxh = 0.04
    maxh_edge = 0.01
    pnts = [(0,0),(0.6,0,0.01),(0.6,0.2,0.003),(3,0.2),(3,1),(0,1)]
else:
    maxh = 0.1
    maxh_edge = maxh
    pnts = [(0,0),(0.6,0,0.03),(0.6,0.2,0.003),(3,0.2),(3,1),(0,1)]
bcs = ["inflow","reflect","reflect","reflect","outflow","reflect"]

pind = [geom.AddPoint(*pnt) for pnt in pnts]
for i in range(len(pind)):
    if(i==2):
        geom.Append(['line',pind[i-1],pind[i]],bc=bcs[i],maxh=maxh_edge)
    else:
        geom.Append(['line',pind[i-1],pind[i]],bc=bcs[i])

mesh = Mesh(geom.GenerateMesh(maxh=maxh,grading=0.15))

dt = 0.025
tend = 1

# using causality constant
local_ctau = True
global_ctau = 2/3
wavespeed = 6
ts = TentSlab(mesh, method="edge")
ts.SetMaxWavespeed(wavespeed)
ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())

order = 4
V = L2(mesh, order=order, dim=mesh.dim+2)
u = GridFunction(V,"u")
cl = Euler(u, ts,
           inflow=mesh.Boundaries("inflow"),
           outflow=mesh.Boundaries("outflow"),
           reflect=mesh.Boundaries("reflect"))
cl.SetTentSolver("SARK",substeps=order*order)

d = 5
rho0 = CoefficientFunction(1.4)
u0 = CoefficientFunction((3,0))
p0 = CoefficientFunction(1)
T0 = 2*p0/rho0
E0 = rho0*(d/4*T0 + 1/2*(u0*u0))

cf = CoefficientFunction((rho0,rho0*u0,E0))
cl.SetInitial(cf)

Draw(u)
visoptions.scalfunction = "u:1"

t = 0
cnt = 0
import time
t1 = time.time()
with TaskManager():
    while t < tend-dt/2:
        cl.Propagate()
        t += dt
        cnt += 1
        if cnt%1 == 0:
            print("{:5f}".format(t))
            Redraw(True)
print(time.time()-t1)
