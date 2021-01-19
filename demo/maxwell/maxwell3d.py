from ngstents.conslaw import Maxwell
from ngstents import TentSlab
from ngsolve.internal import visoptions, viewoptions
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from ngsolve import (Mesh, Draw, Redraw, CoefficientFunction, sqrt, sin, cos,
                     x, y, z, TaskManager)

geom = CSGeometry()
brick = OrthoBrick(Pnt(0, 0, 0), Pnt(1, 1, 1)).bc(2)
geom.Add(brick)
mesh = geom.GenerateMesh(maxh=0.1)
mesh = Mesh(mesh)
# mesh = Mesh("maxwell_cube.vol.gz")

# setting the problem
t_end = 0.5
dt = 0.1
wavespeed = 1

# using causality constant
local_ctau = True
global_ctau = 1/2
ts = TentSlab(mesh, method="edge")#, heapsize=heapsize)
ts.SetWavespeed(wavespeed)
ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())

# setting the approximation space order
order = 2
cl = Maxwell(ts,order)

# initial data
Ez = exp(-40*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5)))
E = CoefficientFunction((0,0,Ez))
H = CoefficientFunction((0,0,0))

cf = CoefficientFunction((E,H))
cl.SetInitial(cf)

sol = cl.sol
cfE = CoefficientFunction((sol[0],sol[1],sol[2]))
cfH = CoefficientFunction((sol[3],sol[4],sol[5]))

Draw(cfE, mesh, name="E")
Draw(cfH, mesh, name="H")

# Draw(sol,sd=1)
visoptions.scalfunction = "E:3"
# visoptions.scalfunction = "u:2"
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
        cl.PropagateSAT(sol.vec, stages=order+1, substeps=2*order)
        t += dt
        cnt += 1
        if cnt%redraw == 0:
            print("{:5f}".format(t))
            # cl.UpdateVisual(sol.vec)
        Redraw(True)
print("total time = ",time.time()-t1)
