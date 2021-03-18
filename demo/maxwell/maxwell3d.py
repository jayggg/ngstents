from ngstents.conslaw import Maxwell
from ngstents import TentSlab
from ngsolve.internal import visoptions, viewoptions
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from ngsolve import (Mesh, Draw, Redraw, CoefficientFunction, sqrt, sin, cos,
                     exp, x, y, z, TaskManager, L2, GridFunction)

geom = CSGeometry()
brick = OrthoBrick(Pnt(0, 0, 0), Pnt(1, 1, 1)).bc("reflect")
geom.Add(brick)
mesh = geom.GenerateMesh(maxh=0.1)
mesh = Mesh(mesh)

# setting the problem
t_end = 0.5
dt = 0.1
wavespeed = 1

# using causality constant
local_ctau = True
global_ctau = 1/2
ts = TentSlab(mesh, method="edge")
ts.SetMaxWavespeed(wavespeed)
ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())

# setting the approximation space order
order = 2
V = L2(mesh, order=order, dim=2*mesh.dim)
u = GridFunction(V,"u")
cl = Maxwell(u, ts, reflect=mesh.Boundaries("reflect"))
cl.SetTentSolver("SAT",stages=order+1, substeps=2*order)

# initial data
Ez = exp(-40*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5)))
E = CoefficientFunction((0,0,Ez))
H = CoefficientFunction((0,0,0))

cf = CoefficientFunction((E,H))
cl.SetInitial(cf)

cfE = CoefficientFunction((u[0],u[1],u[2]))
cfH = CoefficientFunction((u[3],u[4],u[5]))

Draw(cfE, mesh, name="E")
Draw(cfH, mesh, name="H")

visoptions.scalfunction = "E:3"
viewoptions.clipping.enable = 1
visoptions.clipsolution = 'scal'
viewoptions.clipping.ny = 0
viewoptions.clipping.nz = -1

redraw = 1
t = 0
cnt = 0
import time
t1 = time.time()
with TaskManager():
    while t < t_end-dt/2:
        cl.Propagate()
        t += dt
        cnt += 1
        if cnt%redraw == 0:
            print("{:5f}".format(t))
            Redraw(True)
print("total time = ",time.time()-t1)
