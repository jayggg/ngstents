from netgen.geom2d import SplineGeometry
from ngsolve import Mesh, Draw, Redraw
from ngsolve import (CoefficientFunction, IfPos, cos, x, y, InnerProduct,
                     OuterProduct, Id)
from ngsolve import L2, GridFunction, TaskManager
from ngsolve import specialcf as scf
from ngsolve.internal import visoptions
from ngstents import TentSlab
from ngstents.conslaw import ConservationLaw
from math import pi
import time

# define geometry
a = 0.25
geom = SplineGeometry()
pnts = [(-1, 0), (-a, 0), (-a, a), (0, a), (a, a), (a, 0), (1, 0), (1, 1),
        (-1, 1)]
segs = [(0, 1), (1, 5), (5, 6), (6, 7), (7, 8), (8, 0), (1, 2, 3), (3, 4, 5)]
bcs = ["bottom", "bottom", "bottom", "right", "top", "left", "cyl", "cyl"]
dominout = [None, (2, 0), None, None, None, None, (1, 2), (1, 2)]
pind = [geom.AddPoint(*pnt) for pnt in pnts]
for i, seg in enumerate(segs):
    if (dominout[i]):
        leftdomain = dominout[i][0]
        rightdomain = dominout[i][1]
    else:
        leftdomain = 1
        rightdomain = 0
    if (len(seg) == 3):
        geom.Append(['spline3', *seg],
                    leftdomain=leftdomain,
                    rightdomain=rightdomain)
    else:
        geom.Append(['line', *seg],
                    bc=bcs[i],
                    leftdomain=leftdomain,
                    rightdomain=rightdomain)
maxh = 0.1
geom.SetDomainMaxH(2, maxh / 2)
ngmesh = geom.GenerateMesh(maxh=maxh)
mesh = Mesh(ngmesh)
mesh.Curve(3)
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

tend = 2.5
dt = 0.025
# define wave speed
# c = 1   for |x| > a (domain1)
# c = 1/2 for |x| < a (domain2)
wavespeed = CoefficientFunction([1, 1 / 2])

# using causality constant
local_ctau = True
global_ctau = 1 / 4
ts = TentSlab(mesh, method="edge")
ts.SetMaxWavespeed(wavespeed)
ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())

# vector-valued coefficient function
n = scf.normal(mesh.dim)


def Flux(u):
    """
    Compute the flux f(u) for given TrialFunction u,
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
    solves "y = g(u) - (f(u),gradphi)" for u
    g(u) = [[1/c^2 0],[0, 1]] * u
    """
    norm_sqr = InnerProduct(ts.gradphi, ts.gradphi)
    ip = InnerProduct(y[0:mesh.dim], ts.gradphi)
    mu = (y[mesh.dim] + wavespeed**2 * ip) / (1 - wavespeed**2 * norm_sqr)
    q = wavespeed**2 * (y[0:mesh.dim] + mu * ts.gradphi)
    return CoefficientFunction((q, mu))


def ReflectBnd(u):
    """
    reflects values of u at the boundary using the normal vector n
    """
    q = u[0:mesh.dim]
    mu = u[mesh.dim]
    return CoefficientFunction((q - 2 * OuterProduct(n, n) * q, mu))


def BndNumFlux(um):
    """
    defines numerical flux on boundary elements using the normal vector n
    um: trace of u for current element on facet
    """
    # set (q,n) = 0
    return CoefficientFunction((um[mesh.dim] * n, 0))


# define conservation law
order = 2
V = L2(mesh, order=order, dim=mesh.dim + 1)
gfu = GridFunction(V, name="u")
cl = ConservationLaw(gfu,
                     ts,
                     flux=Flux,
                     numflux=NumFlux,
                     inversemap=InverseMap)
cl.SetTentSolver("SARK", substeps=4 * order)
# cl.SetTentSolver("SAT", stages=order+1, substeps=4*order)

# set inital data
mu0 = CoefficientFunction(0)
q0 = CoefficientFunction((0, 0))
cl.SetInitial(CoefficientFunction((q0, mu0)))


# define plane wave for boundary condition
def f(s):
    return (CoefficientFunction(1) - cos(4 * pi * (s - 1))) * IfPos(
        s - 1, IfPos(3 / 2 - s, 1, 0), 0)


k = CoefficientFunction((1, 0))
pos = CoefficientFunction((x, y))


def uex(time):
    return CoefficientFunction((k, 1)) * f(time - InnerProduct(k, pos))


cl.SetBoundaryCF(
    mesh.BoundaryCF({
        "left": NumFlux(cl.u_minus, uex(cl.time)),
        "bottom|right|top": BndNumFlux(cl.u_minus)
    }))

Draw(gfu)
visoptions.scalfunction = "u:3"

t = 0
cnt = 0
redraw = 1

input("press enter to start")
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
