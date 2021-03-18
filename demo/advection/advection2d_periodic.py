from ngsolve import Mesh, Draw, Redraw
from ngsolve import CoefficientFunction, exp, sqrt, sin, cos, x, y, z
from ngsolve import L2, GridFunction
from ngsolve import TaskManager
from ngsolve.internal import visoptions
from ngstents import TentSlab
from ngstents.conslaw import Advection
import time

def Make2DPeriodicMesh(xint, yint, maxh):
    # generate periodic mesh for xint x yint = [xmin, xmax] x [ymin, ymax]
    from netgen.geom2d import SplineGeometry
    periodic = SplineGeometry()
    pnts = [(xint[0],yint[0]),(xint[1],yint[0]),(xint[1],yint[1]),(xint[0],yint[1])]
    pnums = [periodic.AppendPoint(*p) for p in pnts]
    lbot = periodic.Append ( ["line", pnums[0], pnums[1]], bc="bottom")
    lright = periodic.Append ( ["line", pnums[1], pnums[2]], bc="right")
    periodic.Append ( ["line", pnums[0], pnums[3]], leftdomain=0, rightdomain=1, bc="left", copy=lright)
    periodic.Append ( ["line", pnums[3], pnums[2]], leftdomain=0, rightdomain=1, bc="top", copy=lbot)
    return periodic.GenerateMesh(maxh=maxh)

maxh = 0.1
mesh = Mesh(Make2DPeriodicMesh([0,1], [0,1], maxh))

# setting the problem
tend = 10
dt = 0.2
wavespeed = 2

# using causality constant
local_ctau = True
global_ctau = 1
ts = TentSlab(mesh, method="edge")
ts.SetMaxWavespeed(wavespeed)
ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())

order = 4
V = L2(mesh, order=order)
u = GridFunction(V)
cl = Advection(u, ts)
flux = (1,0.1)
cl.SetVectorField( CoefficientFunction(flux) )
cl.SetTentSolver("SAT",stages=order+1, substeps=2*order)

pos = (0.5,0.5)
u0 = CoefficientFunction( exp(-100* ((x-pos[0])*(x-pos[0])+(y-pos[1])*(y-pos[1])) ))
cl.SetInitial(u0)

Draw(u)
visoptions.autoscale = 0

t = 0
cnt = 0
redraw = 1

import time
input("press enter to start")
t1 = time.time()
with TaskManager():
    while t < tend-dt/2:
        cl.Propagate()
        t += dt
        cnt += 1
        if cnt%redraw == 0:
            print("{:5f}".format(t))
            Redraw(True)
print("total time = ",time.time()-t1)
