"""
OBSOLETE / NEEDS REVISION ...
"""
from ngsolve import Mesh, Draw, Redraw
from ngsolve import CoefficientFunction, cos, x, InnerProduct, Id
from ngsolve import L2, GridFunction, TaskManager
from ngsolve import specialcf as scf
from ngsolve.internal import visoptions, viewoptions
from ngstents import TentSlab
from ngstents.meshes import Make1DMesh
from ngstents.conslaw import ConservationLaw
from math import pi

N = 128
mesh = Mesh(Make1DMesh([[0, 2], [2, 4]], [N, N], bcname=["left", "right"]))
'''
solve wave equation: d_tt(phi) - div(alpha grad(phi)) = 0
define first order system

d_t(g(u)) - div f(u) = 0

for u = [q,mu]^t with q = -alpha grad(phi) and mu = d_t(phi)
g(u) = [[alpha^-1 0],[0, 1]] * u
f(u) = [[ I*mu ],
        [ q^t  ]]
isentropic material: alpha = I*c^2 for wave speed c
'''

# wave speeds in domains
ws = CoefficientFunction([1, 1 / 2])

tend = 7
dt = 0.1
local_ctau = True
global_ctau = 1 / 2
ts = TentSlab(mesh, method="edge")
ts.SetMaxWavespeed(ws)
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
    """
    norm_sqr = InnerProduct(ts.gradphi, ts.gradphi)
    mu = (y[mesh.dim] + ws**2 * InnerProduct(y[0:mesh.dim], ts.gradphi)) / (
        1 - ws**2 * norm_sqr)
    q = ws**2 * (y[0:mesh.dim] + mu * ts.gradphi)
    return CoefficientFunction((q, mu))


def BndNumFlux_transparent(um):
    """
    defines numerical flux on boundary elements using the normal vector n
    um: trace of u for current element on facet
    """
    return CoefficientFunction(
        (0.5 * (um[mesh.dim] + 1 / ws * um[0:mesh.dim] * n) * n,
         0.5 * (um[mesh.dim] + (2 - 1 / ws) * um[0:mesh.dim] * n)))


cl = ConservationLaw(gfu,
                     ts,
                     flux=Flux,
                     numflux=NumFlux,
                     inversemap=InverseMap)
cl.SetTentSolver("SAT", stages=order + 1, substeps=4 * order)
# cl.SetTentSolver("SARK", substeps=4*order)

# set inital data
cl.SetInitial(CoefficientFunction((0, 0)))


# define plane wave for boundary condition
def f(s):
    return (CoefficientFunction(1) - cos(4 * pi * (s))) * IfPos(
        s, IfPos(1 / 2 - s, 1, 0), 0)


def planewave(t):
    return CoefficientFunction((ws**2, ws)) * f(t - x)


t = cl.time  # advancing front
cl.SetBoundaryCF(
    mesh.BoundaryCF({
        "left": NumFlux(cl.u_minus, planewave(t)),
        "right": BndNumFlux_transparent(cl.u_minus)
    }))

Draw(gfu)
visoptions.scalfunction = "u:2"
viewoptions.drawedges = 1

redraw = 1
t = 0
cnt = 0
input("start")
import time

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
