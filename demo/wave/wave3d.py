import time
from math import pi
from ngstents.conslaw import Wave
from ngstents import TentSlab
from ngsolve.internal import visoptions, viewoptions
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from ngsolve import (Mesh, Draw, Redraw, CoefficientFunction, sqrt, sin, cos,
                     x, y, z, TaskManager, L2, GridFunction, Integrate, InnerProduct)


geom = CSGeometry()
brick = OrthoBrick(Pnt(0, 0, 0), Pnt(pi, pi, pi)).bc("reflect")
geom.Add(brick)
mesh = geom.GenerateMesh(maxh=0.5)
mesh = Mesh(mesh)

# setting the problem
period = sqrt(3)*pi
n_periods = 1
t_end = period*n_periods
dt = t_end/12
heapsize = 10*1000*1000
wavespeed = 1

# using causality constant
local_ctau = True
global_ctau = 1/2
ts = TentSlab(mesh, method="edge", heapsize=heapsize)
ts.SetMaxWavespeed(wavespeed)
ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())

# setting the approximation space order
order = 2
V = L2(mesh, order=order, dim=mesh.dim+1)
u = GridFunction(V,"u")
wave = Wave(u, ts, reflect=mesh.Boundaries("reflect"))
wave.SetTentSolver("SAT", stages=order+1, substeps=4*order)

mu0 = CoefficientFunction(cos(x)*cos(y)*cos(z))
q0 = CoefficientFunction((0, 0, 0))
cf = CoefficientFunction((q0, mu0))
wave.SetInitial(cf)

Draw(u)
visoptions.scalfunction = "u:4"
viewoptions.clipping.enable = 1
visoptions.clipsolution = 'scal'
viewoptions.clipping.dist = 0.5

t = 0
cnt = 0

redraw = 1
input('start')
t1 = time.time()
with TaskManager():
    print("starting...")
    while t < t_end - dt/2:
        wave.Propagate()
        t += dt
        cnt += 1
        if cnt % redraw == 0:
            print("{:.3f}".format(t))
            Redraw(True)
print("total time = {}".format(time.time()-t1))
# exact solution
exsol = CoefficientFunction((sin(x)*cos(y)*cos(z)*sin(sqrt(3)*t_end)/sqrt(3),
                             cos(x)*sin(y)*cos(z)*sin(sqrt(3)*t_end)/sqrt(3),
                             cos(x)*cos(y)*sin(z)*sin(sqrt(3)*t_end)/sqrt(3),
                             cos(x)*cos(y)*cos(z)*cos(sqrt(3)*t_end)))
Draw(exsol, mesh, 'exact')
l2error = sqrt(Integrate(InnerProduct(u-exsol, u-exsol), mesh, order=3*order))
print("l2error = ", l2error)
