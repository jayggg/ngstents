from netgen.geom2d import SplineGeometry
from ngsolve import Mesh, Draw, Redraw
from ngsolve import CoefficientFunction, sqrt, sin, cos, x, y, z
from ngsolve import TaskManager
from ngsolve.internal import visoptions, viewoptions
from ngstents import TentSlab
from ngstents.conslaw import Wave
from math import pi
import time


geom = SplineGeometry()
geom.AddRectangle(p1=(0,0),p2=(pi,pi),bc=2)
mesh = geom.GenerateMesh(maxh=0.25)
mesh = Mesh(mesh)


# setting the problem
period = sqrt(2)*pi
n_periods = 1
t_end = period*n_periods
dt = t_end/10
heapsize = 10*1000*1000
wavespeed = 1

# using causality constant
local_ctau = True
global_ctau = 2/3
ts = TentSlab(mesh, method="edge", heapsize=heapsize)
ts.SetWavespeed(wavespeed)
ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())

# setting the approximation space order
order = 2
wave = Wave(ts, order)

p = CoefficientFunction(cos(x)*cos(y))
u = CoefficientFunction((0,0))
cf = CoefficientFunction((u, p))
wave.SetInitial(cf)
sol = wave.sol
Draw(sol)
visoptions.scalfunction = "u:4"

t = 0
cnt = 0

redraw = 1
input('start')
t1 = time.time()
with TaskManager():
    print("starting...")
    while t < t_end - dt/2:
        wave.PropagateSAT(sol.vec, stages=order+1, substeps=4*order)
        t += dt
        cnt += 1
        if cnt % redraw == 0:
            print("{:.3f}".format(t))
            Redraw(True)
print("total time = {}".format(time.time()-t1))

# exact solution
exsol = CoefficientFunction((sin(x)*cos(y)*sin(sqrt(2)*t_end)/sqrt(2),
                             cos(x)*sin(y)*sin(sqrt(2)*t_end)/sqrt(2),
                             cos(x)*cos(y)*cos(sqrt(2)*t_end)))
Draw(exsol, mesh, 'exact')
l2error = sqrt(Integrate(InnerProduct(sol-exsol,sol-exsol),mesh,order=3*order))
print("l2error = ",l2error)
