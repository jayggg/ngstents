"""
Burgers equation in 2-D
"""
from netgen.geom2d import SplineGeometry
from ngsolve import Mesh, CoefficientFunction, x, y, exp, Draw, Redraw
from ngsolve import TaskManager, SetNumThreads
from tents import ConsLaw

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

order = 4

cl = ConsLaw(mesh, "burgers", order=order)

sol = cl.sol

cf = CoefficientFunction(exp(-50*((x-0.3)*(x-0.3)+(y-0.3)*(y-0.3))))

cl.SetInitial(cf)

Draw(sol)  # ,sd=5,autoscale=False)

dt = 0.05
tend = 3*dt
t = 0
cnt = 0
cl.CreateTents(dt, 16)
print("max slope", cl.MaxSlope())
if ngs_gui:
    Draw(cl)

if vtk_tents:
    cl.VTKTents('temp')

input('start')
with TaskManager():
    while t < tend - dt/2:
        cl.PropagatePicard(sol.vec, steps=order*order)
        t += dt
        cnt += 1
        if cnt % 1 == 0:
            print("{:.3f}".format(t))
            Redraw(True)
            if step_sol:
                input('step')

if sol_vec:
    print(sol.vec)
