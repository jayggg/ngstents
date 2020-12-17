"""
Burgers equation in 2-D
"""
from netgen.geom2d import SplineGeometry
from ngsolve import Mesh, CoefficientFunction, x, y, exp, Draw, Redraw
from ngsolve import TaskManager, SetNumThreads
from ngstents import TentSlab
from ngstents.conslaw import Burgers

# options
saved_mesh = False
step_sol = False
sol_vec = False
vtk_tents = False
ngs_gui = False

SetNumThreads(1)
if saved_mesh:
    mesh = Mesh("../mesh_geom/square_outflow.vol.gz")
else:
    maxh = 0.05
    geom = SplineGeometry()
    geom.AddRectangle((0, 0), (1, 1), bc=1)
    mesh = Mesh(geom.GenerateMesh(maxh=maxh))

heapsize = 10*1000*1000
dt = 0.05
c = 6
ts = TentSlab(mesh, heapsize=heapsize)
ts.SetWavespeed(c)
ts.PitchTents(dt)
print("max slope", ts.MaxSlope())

order = 4
burg = Burgers(ts, order=order)

sol = burg.sol

cf = CoefficientFunction(exp(-50*((x-0.3)*(x-0.3)+(y-0.3)*(y-0.3))))

burg.SetInitial(cf)

Draw(sol)  # ,sd=5,autoscale=False)

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
        burg.PropagateSARK(sol.vec, substeps=order*order)
        t += dt
        cnt += 1
        if cnt % 1 == 0:
            print("{:.3f}".format(t))
            Redraw(True)
            if step_sol:
                input('step')

if sol_vec:
    print(sol.vec)
