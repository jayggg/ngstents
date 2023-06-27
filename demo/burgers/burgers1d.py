"""
Solve the Burgers equation in 1-D
(Use 1D Netgen Gui:  netgen burgers1d.py)
"""

from ngstents import TentSlab
from ngstents.conslaw import Burgers
from ngsolve.meshes import Make1DMesh
from ngsolve import x, exp, Draw, Redraw
from ngsolve import L2, GridFunction, TaskManager
from ngsolve.internal import viewoptions

step_sol = False  # Step through the solution

dt = 0.05  # tent slab height
c = 4.0  # maximal wave speed
p = 4  # spatial polynomial degree of approximation

mesh = Make1DMesh(100)

ts = TentSlab(mesh)
ts.SetMaxWavespeed(c)
ts.PitchTents(dt)
print("Max Slope: ", ts.MaxSlope())
#  ts.DrawPitchedTentsPlt()

V = L2(mesh, order=p)
u = GridFunction(V, "u")
burg = Burgers(u,
               ts,
               inflow=mesh.Boundaries("left"),
               outflow=mesh.Boundaries("right"))
burg.SetTentSolver(substeps=p * p)
burg.SetInitial(0.5 * exp(-100 * (x - 0.2) * (x - 0.2)))
Draw(u)

tend = 4
t = 0
cnt = 0

viewoptions.drawedges = 1
viewoptions.drawcolorbar = 0

input('start')
with TaskManager():
    while t < tend:
        burg.Propagate()
        t += dt
        cnt += 1
        if cnt % 5 == 0:
            print("{0:.5f}".format(t))
            Redraw()
            if step_sol:
                input('step')
