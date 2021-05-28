from netgen.geom2d import SplineGeometry, unit_square
from ngsolve import Mesh, Draw, Redraw
from ngsolve import (CoefficientFunction, IfPos, sqrt, sin, cos, exp, log, x, y, z,
                     InnerProduct, Norm, OuterProduct, Id)
from ngsolve import L2, GridFunction, TaskManager, Timers, SetNumThreads, Integrate
from ngsolve import specialcf as scf
from ngsolve.internal import visoptions, viewoptions
from ngstents import TentSlab
from ngstents.conslaw import ConservationLaw
from ngstents.utils import Make1DMesh
from math import pi
import time

# define geometry
maxh = 0.05
mesh = Mesh(unit_square.GenerateMesh(maxh=maxh))

'''
solve euler equation:

d_t(g(u)) - div f(u) = 0

for u = [rho, m, E]^t
g(u) = u
f(u) = [[ m^t ],
        [ p*I + m^t*m/rho ],
        [ (E + p)*m^t/rho]]

p = 1/2*rho*T
T = 4/d (E/rho - 1/2 ||m||^2/rho^2)
degrees of freedom of the gas particle d = 5
'''

tend = 0.25
dt = 0.05

# using causality constant
local_ctau = True
global_ctau = 2/3
ts = TentSlab(mesh, method="edge")
ts.SetMaxWavespeed(6)
ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())

# vector-valued coefficient function
n = scf.normal(mesh.dim)
h = scf.mesh_size

def Flux(u):
    """
    Compute the flux f(u) for given TrialFunction u,
    a matrix-valued CoefficientFunction of size (h,w) where
    h represents the number of equations and
    w should correspond to the dimension of the mesh/space
    """
    m = u[1:mesh.dim+1]
    p = 2/5 * (u[mesh.dim+1] - 1/2 * InnerProduct(m,m)/u[0])
    return CoefficientFunction((m,
                                Id(mesh.dim)*p + OuterProduct(m,m)/u[0],
                                (u[mesh.dim+1] + p)/u[0] * m),
                               dims=(V.dim, mesh.dim))

def NumFlux(um, up):
    """
    Compute the numeric flux
    Inputs
    um: trace of u for current element on facet
    up: trace of u for other element on facet
    """
    # um.dims = mesh.dim+2
    # upwind flux
    flux = 0.5 * (Flux(um) + Flux(up))*n + (um - up)
    return flux

def InverseMap(y):
    """
    solves "y = u - (f(u),gradphi)" for u
    
    [derived in disseration "Mapped Tent Pitching Schemes
    for Hyperbolic Systems" by C. Wintersteiger (2020)]
    """
    y_rho = y[0]
    y_m = y[1:mesh.dim+1]
    y_E = y[mesh.dim+1]
    d = 5
    a1 = d/2 * ( y_rho - InnerProduct(y_m, ts.gradphi) )
    a2 = 2*y_E*y_rho - InnerProduct(y_m,y_m)
    normgrad = InnerProduct(ts.gradphi,ts.gradphi)
    p = a2 / (a1 + sqrt( a1**2 - (d+1)*normgrad*a2))
    rho = y_rho**2 / (y_rho - (InnerProduct(y_m, ts.gradphi) + p*normgrad))
    m = rho/y_rho * (y_m + p*ts.gradphi)
    E = 1/y_rho * (rho*y_E + p*InnerProduct(m,ts.gradphi))
    return CoefficientFunction((rho, m, E))

def BndNumFlux(um):
    """
    defines numerical flux on boundary elements using the normal vector n
    um: trace of u for current element on facet
    """
    return CoefficientFunction(( um[mesh.dim]*n, 0))

def Entropy(u):
    """
    Compute the entropy E for given TrialFunction u
    """
    d = 5
    rho = u[0]
    m = u[1:mesh.dim+1]
    T = 4/d * (u[mesh.dim+1]/rho - 1/2 * InnerProduct(m,m)/rho**2)
    T = IfPos(-T, T, 1e-10)
    return rho * ( log(rho) - d/2 * log(T) )

def EntropyFlux(u):
    """
    Compute the entropy flux F for given TrialFunction u
    """
    return u[1:mesh.dim+1]/u[0] * Entropy(u)

def NumEntropyFlux(um, up):
    """
    Compute the numerical flux of the entropy flux F for given
    um: trace of u for current element on facet
    up: trace of u for other element on facet
    """
    mn = InnerProduct(um[1:mesh.dim+1],n)
    return IfPos(mn, mn/um[0] * Entropy(um), InnerProduct(up[1:mesh.dim+1],n)/up[0] * Entropy(up))

def ViscosityCoefficient(u, res):
    """
    Compute the viscosity coefficient for given TrialFunction u
    and entropy residual res
    """
    nu_entr = (h/order)**2 * IfPos(res, res, -res)
    d = 5
    rho = u[0]
    ip_m = InnerProduct(u[1:mesh.dim+1],u[1:mesh.dim+1])
    T = 4/d * (u[mesh.dim+1]/rho - 1/2 * ip_m/rho**2)
    nu_max = 1/20 * h/order * ( sqrt(ip_m) + rho * sqrt((d+2)/d * T) )
    return IfPos( nu_max - nu_entr, nu_entr, nu_max)

# define conservation law
order = 2
V = L2(mesh, order=order, dim=mesh.dim+2)
gfu = GridFunction(V,name="u")
cl = ConservationLaw(gfu, ts,
                     flux=Flux, numflux=NumFlux, inversemap=InverseMap,
                     entropy = Entropy, entropyflux = EntropyFlux,
                     numentropyflux = NumEntropyFlux, visccoeff = ViscosityCoefficient,
                     compile = True)
cl.SetTentSolver("SARK", substeps=2*order)

# set inital data
d = 5
rho = CoefficientFunction(0.1+exp(-200*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))))
m = CoefficientFunction((0,0))
p = CoefficientFunction(0.1+exp(-200*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))))
T = 2*p/rho
E = d/4*T*rho + 1/(2*rho)*m*m

cf = CoefficientFunction((rho,m,E))
cl.SetInitial(CoefficientFunction(cf))

cl.SetBoundaryCF(mesh.BoundaryCF({"left|bottom|right|top" : NumFlux(cl.u_minus, cl.u_minus)}))

Draw(gfu)
visoptions.scalfunction = "u:1"

t = 0
cnt = 0
redraw = 1

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
