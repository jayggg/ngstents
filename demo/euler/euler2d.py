from netgen.geom2d import SplineGeometry
from ngsolve import Mesh, Draw, Redraw, CoefficientFunction, sqrt, x, y, TaskManager
from ngsolve.internal import visoptions
from ngstents import TentSlab
from ngstents.conslaw import Euler

maxh = 0.05
geom = SplineGeometry()
geom.AddRectangle((0,0),(1,1), bc=2)
mesh = Mesh(geom.GenerateMesh(maxh=maxh))

dt = 0.05
tend = 0.25

# using causality constant
local_ctau = True
global_ctau = 2/3
wavespeed = 6
ts = TentSlab(mesh, method="edge")
ts.SetWavespeed(wavespeed)
ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())

order = 4
V = L2(mesh, order=order, dim=mesh.dim+2)
u = GridFunction(V,"u")
cl = Euler(u, ts)
cl.SetTentSolver("SARK",substeps=2*order)

d = 5
rho = CoefficientFunction(0.1+exp(-200*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))))
m = CoefficientFunction((0,0))
p = CoefficientFunction(0.1+exp(-200*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))))
T = 2*p/rho
E = d/4*T*rho + 1/(2*rho)*m*m

cf = CoefficientFunction((rho,m,E))
cl.SetInitial(cf)
Draw(u)

visoptions.scalfunction = "u:1"
visoptions.vecfunction = None

t = 0
cnt = 0
with TaskManager():
    while t<tend-dt/2:
        # cl.PropagateSARK(substeps=2*order)
        cl.Propagate()
        t += dt
        cnt += 1
        if cnt%1 == 0:
            print("{0:.5f}".format(t))
            Redraw(True)
