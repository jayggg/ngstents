from netgen.geom2d import SplineGeometry
from ngsolve import Mesh, Draw, Redraw
from ngsolve import CoefficientFunction, exp, sqrt, sin, cos, x, y, z
from ngsolve import TaskManager
from ngsolve.internal import visoptions, viewoptions
from ngstents import TentSlab
from ngstents.conslaw import Advection
from math import pi
import time

geom = SplineGeometry()
geom.AddCircle(c=(0,0),r=1,bc=3)
mesh = geom.GenerateMesh(maxh=0.15)
mesh = Mesh(mesh)

# setting the problem
tend = 10
dt = tend/10
wavespeed = 2

# using causality constant
local_ctau = True
global_ctau = 1
ts = TentSlab(mesh, method="edge")
ts.SetWavespeed(wavespeed)
ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())

order = 4
cl = Advection(ts,order)
cl.SetFluxField( CoefficientFunction((y,-x)) )

u = CoefficientFunction( exp(-100* ((x-0.5)*(x-0.5)+y*y)))
cl.SetInitial(u)

Draw(cl.sol)
visoptions.autoscale = 0

t = 0
cnt = 0
redraw = 1

import time
# input("press enter to start")
t1 = time.time()
with TaskManager():
    while t < tend-dt/2:
        cl.PropagateSAT(cl.sol.vec,substeps=2*order,stages=order+1)
        t += dt
        cnt += 1
        if cnt%redraw == 0:
            print("{:5f}".format(t))
        Redraw(True)
print("total time = ",time.time()-t1)

# exact solution
exsol = CoefficientFunction(exp(-100*(
    (x-0.5*cos(tend))
    *(x-0.5*cos(tend))
    +(y+0.5*sin(tend))
    *(y+0.5*sin(tend))
)))

Draw(exsol,mesh,'exact')
error = cl.sol - exsol
Draw(error,mesh,'error')
visoptions.scalfunction = 'u:'

# 
l2error = Integrate(error*error,mesh,order=3*order)
print("l2error = ",sqrt(l2error))
print("nodf = ",cl.space.ndof)
