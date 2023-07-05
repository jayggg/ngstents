from netgen.geom2d import SplineGeometry
from ngsolve import (Mesh, Draw, Redraw, CoefficientFunction, IfPos, sqrt, sin,
                     cos, exp, x, y, TaskManager, L2, GridFunction, Integrate,
                     InnerProduct)
from ngsolve.internal import visoptions
from ngstents import TentSlab
from ngstents.conslaw import Wave
import time

geom = SplineGeometry()
geom.AddRectangle(p1=(0, 0), p2=(2, 2), bc="square")
mesh = geom.GenerateMesh(maxh=0.25)
mesh = Mesh(mesh)

# setting the problem
t_end = 2
dt = t_end / 10
wavespeed = 1

# using causality constant
local_ctau = True
global_ctau = 2 / 3
ts = TentSlab(mesh, method="edge")
ts.SetMaxWavespeed(wavespeed)
ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())

# setting the approximation space order
order = 2
V = L2(mesh, order=order, dim=mesh.dim + 1)
u = GridFunction(V, "u")
wave = Wave(u, ts)
wave.SetTentSolver("SAT", stages=order + 1, substeps=4 * order)
# wave.SetTentSolver("SARK", substeps=4*order)


def f(s):
    return (exp(-10 * (s - 1) * (s - 1)) - exp(-10)) / (1 - exp(-10)) * IfPos(
        s, IfPos(2 - s, 1, 0), 0)


k = CoefficientFunction((cos(1), sin(1)))
pos = CoefficientFunction((x, y))


def uex(time):
    return CoefficientFunction((k, 1)) * f(time - InnerProduct(k, pos))


mu0 = CoefficientFunction(0)
q0 = CoefficientFunction((0, 0))
cf = CoefficientFunction((q0, mu0))
wave.SetInitial(cf)

t = wave.time  # advancing front
wave.SetBoundaryCF(mesh.BoundaryCF({"square": uex(t)}))

Draw(u)
visoptions.scalfunction = "u:3"

t = 0
cnt = 0

redraw = 1
# input("start")
t1 = time.time()
with TaskManager():
    while t < t_end - dt / 2:
        wave.Propagate()
        t += dt
        cnt += 1
        if cnt % redraw == 0:
            print("{:.3f}".format(t))
            Redraw(True)
print("total time = {}".format(time.time() - t1))

# exact solution
uexact = uex(t_end)
Draw(uexact, mesh, "u_exact")
l2error = sqrt(
    Integrate(InnerProduct(u - uexact, u - uexact), mesh, order=3 * order))
print("l2error = ", l2error)
