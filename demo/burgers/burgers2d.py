"""
Burgers equation in 2-D
"""
from netgen.geom2d import unit_square
from ngsolve import Mesh, CoefficientFunction, x, y, exp, Draw, Redraw
from ngsolve import L2, GridFunction
from ngsolve import TaskManager, SetNumThreads
from ngstents import TentSlab
from ngstents.conslaw import Burgers

# options
step_sol = False
sol_vec = False
vtk_tents = False
ngs_gui = False

maxh = 0.05
mesh = Mesh(unit_square.GenerateMesh(maxh=maxh))

heapsize = 10*1000*1000
dt = 0.05
c = 6
ts = TentSlab(mesh, heapsize=heapsize)
ts.SetMaxWavespeed(c)
ts.PitchTents(dt)
print("max slope", ts.MaxSlope())

order = 4
V = L2(mesh, order=order)
u = GridFunction(V,"u")
burg = Burgers(u, ts, outflow=mesh.Boundaries("left|bottom|right|top"))
burg.SetTentSolver("SARK",substeps=order*order)

cf = CoefficientFunction(exp(-50*((x-0.3)*(x-0.3)+(y-0.3)*(y-0.3))))
burg.SetInitial(cf)

Draw(u)

tend = 10*dt
t = 0
cnt = 0
if ngs_gui:
    Draw(burg)

if vtk_tents:
    burg.VTKTents('temp')

input('start')
with TaskManager():
    while t < tend - dt/2:
        burg.Propagate()
        t += dt
        cnt += 1
        if cnt % 1 == 0:
            print("{:.3f}".format(t))
            Redraw(True)
            if step_sol:
                input('step')

if sol_vec:
    print(sol.vec)
