from netgen.geom2d import SplineGeometry
from ngsolve import Mesh, Draw, Redraw
from ngsolve import CoefficientFunction, sqrt, sin, cos, x, y
from ngsolve import InnerProduct, Id, L2, GridFunction, TaskManager, Integrate
from ngsolve import specialcf as scf
from ngsolve.internal import visoptions, viewoptions
from ngstents import TentSlab
from ngstents.conslaw import ConservationLaw
from math import pi
import time

geom = SplineGeometry()
geom.AddRectangle((0, 0), (pi, pi), bc="reflect")
ngmesh = geom.GenerateMesh(maxh=0.25)
mesh = Mesh(ngmesh)

# setting the problem
tend = 2 / sqrt(2) * pi
dt = tend / 10
wavespeed = 1

# using causality constant
local_ctau = True
global_ctau = 2 / 3
ts = TentSlab(mesh, method="edge")
ts.SetMaxWavespeed(wavespeed)
ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())

order = 2
V = L2(mesh, order=order, dim=mesh.dim + 1)
gfu = GridFunction(V, name="u")

# vector-valued coefficient function
n = scf.normal(mesh.dim)


def Flux(u):
    """
    Compute the flux for given TrialFunction u,
    a matrix-valued CoefficientFunction of size (h,w) where
    h represents the number of equations and
    w should correspond to the dimension of the mesh/space
    """
    return CoefficientFunction((Id(mesh.dim) * u[mesh.dim], u[0:mesh.dim]),
                               dims=(V.dim, mesh.dim))


def NumFlux(um, up):
    """
    Compute the numeric flux
    Inputs
    um: trace of u for current element on facet
    up: trace of u for other element on facet
    """
    # upwind flux
    flux = 0.5 * (Flux(um) + Flux(up)) * n
    flux_vec = 0.5 * (um[0:mesh.dim] - up[0:mesh.dim]) * n * n
    flux_scal = 0.5 * (um[mesh.dim] - up[mesh.dim])
    return flux + CoefficientFunction((flux_vec, flux_scal))


def InverseMap(y):
    """
    solves "y = u - (f(u),gradphi)" for u

    assuming wave speed = 1
    """
    norm_sqr = InnerProduct(ts.gradphi, ts.gradphi)
    mu = (y[mesh.dim] + InnerProduct(y[0:mesh.dim], ts.gradphi)) / (1 -
                                                                    norm_sqr)
    q = y[0:mesh.dim] + mu * ts.gradphi
    return CoefficientFunction((q, mu))


def BndNumFlux(um):
    """
    defines numerical flux on boundary elements using the normal vector n
    um: trace of u for current element on facet
    """
    # set (q,n) = 0
    return CoefficientFunction((um[mesh.dim] * n, 0))


cl = ConservationLaw(gfu,
                     ts,
                     flux=Flux,
                     numflux=NumFlux,
                     inversemap=InverseMap)
cl.SetBoundaryCF(mesh.BoundaryCF({"reflect": BndNumFlux(cl.u_minus)}))
# cl.SetTentSolver("SAT",stages=order+1, substeps=4*order)
cl.SetTentSolver("SARK", stages=order + 1, substeps=4 * order)

mu0 = cos(x) * cos(y)
q0 = CoefficientFunction(tuple(mesh.dim * [0]))
cl.SetInitial(CoefficientFunction((q0, mu0)))

Draw(gfu)
visoptions.scalfunction = "u:{:0}".format(mesh.dim + 1)
viewoptions.drawedges = 1

t = 0
cnt = 0
redraw = 1

# input("press enter to start")
t1 = time.time()
with TaskManager():
    while t < tend - dt / 2:
        cl.Propagate()
        t += dt
        cnt += 1
        if cnt % redraw == 0:
            print("{:5f}".format(t))
            Redraw(True)
print("total time = ", time.time() - t1)

exsol = CoefficientFunction((sin(x) * cos(y) * sin(sqrt(2) * tend) / sqrt(2),
                             cos(x) * sin(y) * sin(sqrt(2) * tend) / sqrt(2),
                             cos(x) * cos(y) * cos(sqrt(2) * tend)))

Draw(exsol, mesh, "exact")
l2error = sqrt(
    Integrate(InnerProduct(gfu - exsol, gfu - exsol), mesh, order=3 * order))
print("l2error = ", l2error)
