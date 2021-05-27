"""
Solve the Burgers equation in 1-D

Note: zoom in Netgen to see solution for cfnum = 2
"""
from ngstents import TentSlab
from ngstents.utils import Make1DMesh, Make1DMeshSpecified
from ngstents.conslaw import Burgers
from ngsolve import Mesh, CoefficientFunction, x, exp, Draw, Redraw
from ngsolve import L2, GridFunction
from ngsolve import TaskManager, SetNumThreads, Timers
from ngsolve.internal import viewoptions

SetNumThreads(1)

# Options
cfnum = 0
draw_tents = False  # Plot tents.  Don't set this in Netgen context.
step_sol = False     # Step through the solution
spec = False         # Use a non-uniform mesh for more interesting tent plot
dt = 0.05           # Tent slab height
c = 4.0             # Characteristic speed
order = 4           # Order of L2 space for solution

if spec:
    pts =[0, 0.04, 0.08, 0.12, 0.18, 0.24, 0.30, 0.35, 0.40, 0.45, 0.50,
             0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0,]
    mesh = Make1DMeshSpecified(pts, bcname=["left","right"])
else:
    mesh = Make1DMesh([[0, 1]], [100], bcname=["left","right"])
mesh = Mesh(mesh)

ts = TentSlab(mesh)
ts.SetMaxWavespeed(c)
ts.PitchTents(dt)
print("Max Slope: ", ts.MaxSlope())

if draw_tents:
    ts.DrawPitchedTentsPlt()

V = L2(mesh, order=order)
u = GridFunction(V,"u")
burg = Burgers(u, ts, inflow=mesh.Boundaries("left"), outflow=mesh.Boundaries("right"))
burg.SetTentSolver("SARK",substeps=order*order)

cf = CoefficientFunction(0.5*exp(-100*(x-0.2)*(x-0.2))) if cfnum == 0 \
    else CoefficientFunction(x) if cfnum == 1 \
    else CoefficientFunction(1)

burg.SetInitial(cf)
Draw(u)

tend = 4
t = 0
cnt = 0

viewoptions.drawedges = 1
viewoptions.drawcolorbar = 0

input('start')
with TaskManager():
    while t < tend - dt/2:
        burg.Propagate()
        t += dt
        cnt += 1
        if cnt % 1 == 0:
            print("{0:.5f}".format(t))
            Redraw()
            # print("||res|| = ",Norm(burg.res.vec.FV()))
            if step_sol:
                input('step')

    f = open("timings_burgers.csv","w")
for t in Timers():
    if t["name"] == "SARK::Propagate Tent" or \
       t["name"] == "Propagate" or \
       t["name"] == "CalcFluxTent" or \
       t["name"] == "Cyl2Tent" or \
       t["name"] == "ApplyM1" or \
       t["name"] == "Tent2Cyl" or \
       t["name"] == "CalcViscosityTent" or \
       t["name"] == "calc residual" or \
       t["name"] == "calc nu" or \
       t["name"] == "apply viscosity" or \
       t["name"] ==  "CalcEntropy" or \
       t["name"] ==  "EntropyFlux" or \
       t["name"] ==  "Inverse Map" or \
       t["name"] ==  "Flux" or \
       t["name"] ==  "NumFlux" or \
       t["name"] ==  "EntropyViscCoeff":
        f.write(", ".join([t["name"], str(t["time"]) , str(t["counts"])])+"\n")
        print(t)
f.close()
