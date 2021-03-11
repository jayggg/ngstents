from netgen.geom2d import SplineGeometry
from ngsolve import Mesh, Draw, Redraw
from ngsolve import (CoefficientFunction, IfPos, sqrt, sin, cos, exp, x, y, z,
                     InnerProduct, Norm, OuterProduct, Id)
from ngsolve import L2, GridFunction, TaskManager, SetNumThreads
from ngsolve import specialcf as scf
from ngsolve.internal import visoptions, viewoptions
from ngstents import TentSlab
from ngstents.conslaw import ConservationLaw
from ngstents.utils import Make1DMesh
import time

dim = 2
if(dim == 1):
    N = 100
    ngmesh = Make1DMesh([[0,1]], [N], bc=[1,1])
else:
    geom = SplineGeometry()
    geom.AddRectangle((0,0),(1,1), bcs=[3,3,3,3])
    ngmesh = geom.GenerateMesh(maxh=0.05)
mesh = Mesh(ngmesh)

# setting the problem
tend = 0.4
dt = 0.1
wavespeed = 2

# using causality constant
local_ctau = True
ts = TentSlab(mesh, method="edge")
ts.SetWavespeed(wavespeed)
ts.PitchTents(dt=dt, local_ct=local_ctau)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())

order = 2
V = L2(mesh, order=order, dim=mesh.dim+1)
gfu = GridFunction(V,name="u")

# vector-valued coefficient function
n = scf.normal(mesh.dim)

def Flux(u):
    """
    Compute the flux for given TrialFunction u,
    a matrix-valued CoefficientFunction of size (h,w) where
    h represents the number of equations and
    w should correspond to the dimension of the mesh/space
    """
    q = CoefficientFunction(tuple([u[i] for i in range(0,mesh.dim)]),dims=(mesh.dim,1))
    mu = u[mesh.dim]
    return CoefficientFunction((Id(mesh.dim)*mu,q.trans), dims=(V.dim, mesh.dim))

def NumFlux(um, up):
    """
    Compute the numeric flux
    Inputs
    um: trace of u for current element on facet
    up: trace of u for other element on facet
    """
    # upwind flux
    flux = 0.5 * (Flux(um) + Flux(up))*n
    qm = CoefficientFunction(tuple([um[i] for i in range(0,mesh.dim)]),dims=(mesh.dim,1))
    qp = CoefficientFunction(tuple([up[i] for i in range(0,mesh.dim)]),dims=(mesh.dim,1))
    flux_vec = 0.5 * OuterProduct(n,n) * (qm - qp)
    flux_scal = 0.5 * (um[mesh.dim]-up[mesh.dim])
    return flux + CoefficientFunction((flux_vec,flux_scal))

def InverseMap(y):
    """
    solves "y = u - (f(u),gradphi)" for u
    """
    norm_sqr = InnerProduct(ts.gradphi,ts.gradphi)
    y_q = CoefficientFunction(tuple([y[i] for i in range(0,mesh.dim)]),dims=(mesh.dim,1))
    y_mu = y[mesh.dim]
    ip = InnerProduct(y_q,ts.gradphi)
    mu = (y_mu + InnerProduct(y_q,ts.gradphi))/(1-norm_sqr)
    q = y_q + mu*ts.gradphi
    return CoefficientFunction((q,mu))

cl = ConservationLaw(gfu, ts, flux=Flux, numflux=NumFlux, inversemap=InverseMap, compile=True)
cl.SetTentSolver("SAT",stages=order+1, substeps=4*order)

if(dim == 1):
    x0 = 0.5
    q0 = CoefficientFunction(0)
    mu0 = 0.5*exp(-200* ((x-x0)*(x-x0)))
else:
    x0, y0 = 0.5, 0.5
    q0 = CoefficientFunction( (0,0) )
    mu0 = 0.5*exp(-1e3* ((x-x0)*(x-x0) + (y-y0)*(y-y0)))
cl.SetInitial(CoefficientFunction((q0,mu0)))

Draw(gfu)
visoptions.scalfunction = "u:{:0}".format(mesh.dim+1)
viewoptions.drawedges = 1

t = 0
cnt = 0
redraw = 1

import time
input("press enter to start")
t1 = time.time()
with TaskManager():
    while t < tend-dt/2:
        cl.Propagate()
        t += dt
        cnt += 1
        if cnt%redraw == 0:
            print("{:5f}".format(t))
            Redraw(True)
print("total time = ",time.time()-t1)
