"""
Solve the Burgers equation in 1-D

Notes:
1. zoom in Netgen to see solution for cfnum = 2
2. boundary condition numbers 1 = inflow, 3 = outflow
"""
from ngstents import TentSlab
from ngstents.utils import Make1DMesh, Make1DMeshSpecified
from ngstents.conslaw import Burgers
from ngsolve import Mesh, CoefficientFunction, x, exp, Draw, Redraw
from ngsolve import TaskManager, SetNumThreads
from ngsolve.internal import viewoptions

# Options
cfnum = 0
draw_tents = False  # Plot tents.  Don't set this in Netgen context.
step_sol = True     # Step through the solution
spec = True         # Use a non-uniform mesh for more interesting tent plot
dt = 0.05           # Tent slab height
c = 4.0             # Characteristic speed
order = 4           # Order of L2 space for solution

SetNumThreads(1)
if spec:
    pts =[0, 0.04, 0.08, 0.12, 0.18, 0.24, 0.30, 0.35, 0.40, 0.45, 0.50,
             0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0,]
    mesh = Make1DMeshSpecified(pts, bc=[1, 3])
else:
    mesh = Make1DMesh([[0, 1]], [20], bc=[1, 3])
mesh = Mesh(mesh)

ts = TentSlab(mesh)
ts.SetWavespeed(c)
ts.PitchTents(dt)
print("Max Slope: ", ts.MaxSlope())

if draw_tents:
    ts.DrawPitchedTentsPlt()

burg = Burgers(ts, order)
sol = burg.sol

cf = CoefficientFunction(exp(-50*(x-0.3)*(x-0.3))) if cfnum == 0 \
    else CoefficientFunction(x) if cfnum == 1 \
    else CoefficientFunction(1)

burg.SetInitial(cf)
Draw(sol)

tend = 10
t = 0
cnt = 0

viewoptions.drawedges = 1
viewoptions.drawcolorbar = 0

input('start')
with TaskManager():
    while t < tend:
        burg.PropagateSARK(sol.vec, substeps=order*order)
        t += dt
        cnt += 1
        if cnt % 5 == 0:
            print("{0:.5f}".format(t))
            Redraw()
            if step_sol:
                input('step')
