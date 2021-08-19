from ngsolve.internal import visoptions
from ngstents import TentSlab
from ngstents.trefftz import TWave
from ngsolve.TensorProductTools import *
from ngsolve import *
import time
from math import pi


geom = SplineGeometry()
geom.AddRectangle(p1=(0, 0), p2=(pi, pi), bc="neumann")
mesh = Mesh(geom.GenerateMesh(maxh=0.25))

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
ts.SetMaxWavespeed(wavespeed)
ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)

order=5
mu0 = CoefficientFunction(cos(x)*cos(y))
q0 = CoefficientFunction((0, 0))
cf = CoefficientFunction((q0, mu0))
TT=TWave(order,ts,CoefficientFunction(wavespeed))
TT.SetInitial(cf)
TT.SetBoundaryCF(q0)

t1 = time.time()
t = 0
cnt = 0
redraw = 1
V = L2(mesh,order=order-1)**3
u = GridFunction(V,"u")
Draw(u)
visoptions.scalfunction = "u:3"
with TaskManager():
    while t < t_end - dt/2:
        TT.Propagate()
        t += dt
        # cnt += 1
        # if cnt % redraw == 0:
            # print("{:.3f}".format(t))
            # TT.GetWave(u)
            # Redraw(True)
print("total time = {}".format(time.time()-t1))

# exact solution
exsol = CoefficientFunction((sin(x)*cos(y)*sin(sqrt(2)*t_end)/sqrt(2),
                             cos(x)*sin(y)*sin(sqrt(2)*t_end)/sqrt(2),
                             cos(x)*cos(y)*cos(sqrt(2)*t_end)))
TT.GetWave(u)
l2error = sqrt(Integrate(InnerProduct(u-exsol, u-exsol), mesh, order=3*order))
print("l2error = ", l2error)
