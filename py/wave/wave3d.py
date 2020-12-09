from netgen.csg import *
from math import pi

geom = CSGeometry()
brick=OrthoBrick(Pnt(0,0,0),Pnt(pi,pi,pi)).bc(2)
geom.Add(brick)
mesh = geom.GenerateMesh(maxh=0.25)

from ngsolve import *
from ngsolve.internal import visoptions, viewoptions

from conslaw import *

mesh = Mesh(mesh)

cl = ConsLaw(mesh,"wave",2)
cl.SetupDG()

p = CoefficientFunction(sin(x)*sin(y)*sin(z))
u = CoefficientFunction((0,0,0))

cf = CoefficientFunction((u,p))
cl.SetInitial(cf)

sol = cl.sol
Draw(sol,sd=5)
visoptions.scalfunction = "u:4"
viewoptions.clipping.enable = 1
visoptions.clipsolution = 'scal'

period = sqrt(3)*pi
redraw = 1
nperiods = 1
tend = period*nperiods
dt = period/10

t = 0
cnt = 0
cl.CreateTents(dt,10.0)
# cl.DrawTents()

import time
input("press enter to start")
t1 = time.time()
with TaskManager():
    while t < tend-dt/2:
        cl.PropagateRK(sol.vec)
        t += dt
        cnt += 1
        if cnt%redraw == 0:
            print("{:5f}".format(t))
            cl.UpdateVisual(sol.vec)
        Redraw(True)
print("total time = ",time.time()-t1)

# exact solution
exsol = CoefficientFunction((-cos(x)*sin(y)*sin(z)*sin(sqrt(3)*tend)/sqrt(3),
                             -sin(x)*cos(y)*sin(z)*sin(sqrt(3)*tend)/sqrt(3),
                             -sin(x)*cos(y)*cos(z)*sin(sqrt(3)*tend)/sqrt(3),
                             sin(x)*sin(y)*sin(z)*cos(sqrt(3)*tend)))
Draw(exsol,mesh,'exact')
