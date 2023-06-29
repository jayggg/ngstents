'''
Solve Burgers equation: d_t(u) - div f(u) = 0

f(u) = 1/2 * u^2

'''

from ngsolve import Mesh, Draw, Redraw
from ngsolve import CoefficientFunction, x, cos, sqrt, exp, IfPos, Id, InnerProduct
from ngsolve import L2, GridFunction, TaskManager, SetNumThreads
from ngsolve import specialcf as scf
from ngsolve.internal import visoptions, viewoptions
from ngstents import TentSlab
from ngsolve.meshes import Make1DMesh
from ngstents.conslaw import ConservationLaw
from math import pi

N = 20
mesh = Make1DMesh(N)

tend = 4
dt = 0.05
ts = TentSlab(mesh, method="edge")
ts.SetMaxWavespeed(4)
ts.PitchTents(dt=dt)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())

order = 4
V = L2(mesh, order=order)
gfu = GridFunction(V, name="u")

n = scf.normal(mesh.dim)
h = scf.mesh_size


def Flux(u):
    """
    Compute the flux for given TrialFunction u,
    a matrix-valued CoefficientFunction of size (h,w) where
    h represents the number of equations and
    w should correspond to the dimension of the mesh/space
    """
    return CoefficientFunction(1 / 2 * u**2)


def NumFlux(um, up):
    """
    Compute the numeric flux
    Inputs
    um: trace of u for current element on facet
    up: trace of u for other element on facet
    """
    # upwind flux
    return 0.5 * IfPos((um + up) * n, um**2, up**2) * n


def InverseMap(y):
    """
    solves "y = u - (f(u),gradphi)" for u
    """
    return 2 * y / (1 + sqrt(1 - 2 * ts.gradphi * y))


def Entropy(u):
    """
    Compute the entropy E for given TrialFunction u
    """
    return CoefficientFunction(u**2)


def EntropyFlux(u):
    """
    Compute the entropy flux F for given TrialFunction u
    """
    return CoefficientFunction(1 / 3 * u**3)


def NumEntropyFlux(um, up):
    """
    Compute the numerical flux of the entropy flux F for given
    um: trace of u for current element on facet
    up: trace of u for other element on facet
    """
    umean = 0.5 * (um + up)
    return 1 / 3 * IfPos(umean**2 * n, um**3, up**3) * n


def ViscosityCoefficient(u, res):
    """
    Compute the viscosity coefficient for given TrialFunction u
    and entropy residual res
    """
    nu_entr = (h / order)**2 * IfPos(res, res, -res) / Entropy(u)
    nu_max = h / order * IfPos(u, u, -u)
    return IfPos(nu_max - nu_entr, nu_entr, nu_max)


cl = ConservationLaw(gfu,
                     ts,
                     flux=Flux,
                     numflux=NumFlux,
                     inversemap=InverseMap,
                     entropy=Entropy,
                     entropyflux=EntropyFlux,
                     numentropyflux=NumEntropyFlux,
                     visccoeff=ViscosityCoefficient,
                     compile=True)
cl.SetTentSolver("SARK", stages=3, substeps=order * order)

# set inital data
cf0 = CoefficientFunction(0.5 * exp(-100 * (x - 0.2) * (x - 0.2)))
cl.SetInitial(cf0)

# set boundary conditions
cl.SetBoundaryCF(
    mesh.BoundaryCF({
        "left": NumFlux(cl.u_minus, cf0),
        "right": NumFlux(cl.u_minus, cl.u_minus)
    }))
cl.SetNumEntropyFlux(mesh.BoundaryCF({".*": EntropyFlux(cl.u_minus)}))

Draw(gfu)
viewoptions.drawedges = 1

redraw = 5
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
