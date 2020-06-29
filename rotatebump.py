# PLAN / WISHLIST:
#
# We would like to have tent-based time stepping along
# the lines of this code. Instead of the "ApplyDG" below,
# we would need to make a facility, say "Propagate".
# It should be able to work with numerical/fluxes input
# from python (like below) and be able to perform implicit,
# SAT, and SARK propagations.


from netgen.geom2d import unit_square
from ngsolve import Mesh, CoefficientFunction, exp, x, y
import ngsolve as ng
from tents import ApplyDG

ng.ngsglobals.msg_level = 1

mesh = Mesh(unit_square.GenerateMesh(maxh=0.05))
ng.SetHeapSize(100*1000*1000)
fes = ng.L2(mesh, order=3)
gfu = ng.GridFunction(fes)
gfu.Set(exp(-100 * ((x-0.7)*(x-0.7)+(y-0.5)*(y-0.5))))


b = CoefficientFunction((y-0.5, 0.5-x))
n = ng.specialcf.normal(mesh.dim)


def Flux(u):
    return CoefficientFunction(b*u, dims=(1, 2))


def NumFlux(um, up):
    bn = b*n
    return ng.IfPos(bn, bn*um, bn*up)


ng.Draw(gfu, min=0, max=1)
Au = gfu.vec.CreateVector()

dt = 0.001
tend = 5

t = 0

with ng.TaskManager():
    while t < tend:
        print("t = ", t)
        ApplyDG(fes, Flux, NumFlux, gfu.vec, Au)
        fes.SolveM(rho=CoefficientFunction(1), vec=Au)
        gfu.vec.data += dt * Au
        t += dt
        ng.Redraw(blocking=True)


ng.Draw(Flux(gfu), mesh, "flux")
