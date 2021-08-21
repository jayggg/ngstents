from ngsolve import Mesh, Draw, Redraw
from ngsolve import CoefficientFunction, IfPos, exp, x, y, InnerProduct
from ngsolve import L2, GridFunction, TaskManager
from ngsolve import specialcf as scf
from ngsolve.internal import visoptions
from ngstents import TentSlab
from ngstents.conslaw import ConservationLaw
import time


def Make2DPeriodicMesh(xint, yint, maxh):
    # generate periodic mesh for xint x yint = [xmin, xmax] x [ymin, ymax]
    from netgen.geom2d import SplineGeometry
    periodic = SplineGeometry()
    pnts = [(xint[0], yint[0]), (xint[1], yint[0]),
            (xint[1], yint[1]), (xint[0], yint[1])]
    pnums = [periodic.AppendPoint(*p) for p in pnts]
    lbot = periodic.Append(["line", pnums[0], pnums[1]], bc="bottom")
    lright = periodic.Append(["line", pnums[1], pnums[2]], bc="right")
    periodic.Append(["line", pnums[0], pnums[3]], leftdomain=0,
                    rightdomain=1, bc="left", copy=lright)
    periodic.Append(["line", pnums[3], pnums[2]], leftdomain=0,
                    rightdomain=1, bc="top", copy=lbot)
    return periodic.GenerateMesh(maxh=maxh)


maxh = 0.1
mesh = Mesh(Make2DPeriodicMesh([0, 1], [0, 1], maxh))

# setting the problem
tend = 1
dt = 0.2
wavespeed = 2

# using causality constant
local_ctau = True
global_ctau = 1
ts = TentSlab(mesh, method="edge")
ts.SetMaxWavespeed(wavespeed)
ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())

order = 4
V = L2(mesh, order=order)
u = GridFunction(V)

# vector-valued coefficient function
b = CoefficientFunction((1, 0.1))
n = scf.normal(mesh.dim)


def Flux(u):
    """
    Compute the flux for given TrialFunction u,
    a matrix-valued CoefficientFunction of size (h,w) where
    h represents the number of equations and
    w should correspond to the dimension of the mesh/space
    """
    return CoefficientFunction(b*u, dims=(V.dim, mesh.dim))


def NumFlux(um, up):
    """
    Compute the numeric flux
    Inputs
    um: trace of u for current element on facet
    up: trace of u for other element on facet
    """
    # upwind flux
    bn = b*n
    return IfPos(bn, bn*um, bn*up)


def InverseMap(y):
    """
    solves "y = u - (b*u,gradphi)" for u
    """
    return y/(1-InnerProduct(b, ts.gradphi))


cl = ConservationLaw(u, ts, flux=Flux, numflux=NumFlux, inversemap=InverseMap)

cl.SetTentSolver("SAT", stages=order+1, substeps=2*order)

pos = (0.5, 0.5)
u0 = CoefficientFunction(
    exp(-100 * ((x-pos[0])*(x-pos[0])+(y-pos[1])*(y-pos[1]))))
cl.SetInitial(u0)

Draw(u)
visoptions.autoscale = 0

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
        if cnt % redraw == 0:
            print("{:5f}".format(t))
            Redraw(True)
print("total time = ", time.time()-t1)
