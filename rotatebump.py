# PLAN / WISHLIST:
#
# We would like to have tent-based time stepping along
# the lines of this code. Instead of the "ApplyDG" below,
# we would need to make a facility, say "Propagate".
# It should be able to work with numerical/fluxes input
# from python (like below) and be able to perform implicit,
# SAT, and SARK propagations.

# TODO: Test multiple equations (system)
# TODO: Add global boundary conditions defined in Python
#           inflow, outflow, reflecting, transparent, periodic.
#           User should be able to specify regions
# TODO: Add entropy and entropy flux defined in Python
# TODO: Add computation of viscosity defined in Python
# TODO: Something with tents?  e.g. non-vertical internal facets?
#           supertents?  Not pitching whole slab at once?
# TODO: For propagation, do we want a single "Propagate" facility with
#           options for various methods or would it be possible to
#           specify a sequence of tasks in Python that should be performed
#           each timestep so that variations on the methods could be
#           explored?  Probably don't need this...

from netgen.geom2d import unit_square
from ngsolve import Mesh, CoefficientFunction, exp, x, y
import ngsolve as ng
from tents import ApplyDG

ng.ngsglobals.msg_level = 1
# ng.ngsglobals.pajetrace=10000000

mesh = Mesh(unit_square.GenerateMesh(maxh=0.05))
ng.SetHeapSize(100*1000*1000)

fes = ng.L2(mesh, order=3)
gfu = ng.GridFunction(fes)
gfu.Set(exp(-100 * ((x-0.7)*(x-0.7)+(y-0.5)*(y-0.5))))

# vector-valued coefficient function
b = CoefficientFunction((y-0.5, 0.5-x))
n = ng.specialcf.normal(mesh.dim)


def Flux(u):
    """
    Compute the flux for given TrialFunction u,
    a matrix-valued CoefficientFunction of size (h,w) where
    h represents the number of equations and
    w should correspond to the dimension of the mesh/space
    """
    # print(u) -> coef trial-function diffop = Id, real
    return CoefficientFunction(b*u, dims=(fes.dim, mesh.dim))


def NumFlux(um, up):
    """
    Compute the numeric flux
    Inputs
    um: trace of u for current element on facet
    up: trace of u for other element on facet
    """
    # upwind flux
    bn = b*n
    return ng.IfPos(bn, bn*um, bn*up)


ng.Draw(gfu, min=0, max=1)
Au = gfu.vec.CreateVector()

dt = 0.001

# instability at t = 3.5
tend = 5

t = 0

with ng.TaskManager():
    # Forward Euler timestepping
    while t < tend:
        print("t = ", t)
        # result is returned in Au
        ApplyDG(fes, Flux, NumFlux, gfu.vec, Au)
        # uₖ₊₁ = uₖ + M⁻¹(Au) dt
        fes.SolveM(rho=CoefficientFunction(1), vec=Au)
        gfu.vec.data += dt * Au
        t += dt
        ng.Redraw(blocking=True)


ng.Draw(Flux(gfu), mesh, "flux")
