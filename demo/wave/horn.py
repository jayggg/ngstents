import ngsolve as ng
from netgen.geom2d import EdgeInfo as EI, Solid2d, CSG2d
from ngstents import TentSlab
from ngstents.conslaw import Wave
from ngsolve import L2, GridFunction
from ngsolve.internal import visoptions
from ngsolve import exp, CF, y


geo = CSG2d()

horn = Solid2d(
    [(1, 0.55),
     EI((1, 1), bc='out'),      # right curved boundary (with control point)
     (0, 1),
     EI((-1,  1), bc='out'),    # left curved bdry
     (-1, 0.55),
     EI(bc='cone'),             # conical walls
     (-0.03, 0),
     EI(maxh=0.02, bc='pipe'),  # feed pipe
     (-0.03, -0.5),
     EI(maxh=0.02, bc='in'),
     (+0.03, -0.5),
     EI(maxh=0.02, bc='pipe'),
     (+0.03, 0),
     EI(bc='cone')
     ], mat='air')

geo.Add(horn)
m = geo.GenerateMesh(maxh=0.15)
mesh = ng.Mesh(m)
mesh.Curve(4)
# ng.Draw(mesh)


ts = TentSlab(mesh, method="edge", heapsize=10*1000*1000)
wavespeed = 1
ts.SetMaxWavespeed(wavespeed)
ts.PitchTents(dt=0.1, local_ct=True, global_ct=2/3)
print("max slope", ts.MaxSlope())
print("n tents", ts.GetNTents())


order = 4
V = L2(mesh, order=order, dim=mesh.dim+1)
u = GridFunction(V, "u")
wave = Wave(u, ts,
            inflow=mesh.Boundaries('in'),
            transparent=mesh.Boundaries('out'),
            reflect=mesh.Boundaries('pipe|cone'))

wave.SetTentSolver("SAT", stages=order+1, substeps=4*order)


def f(s):
    d = 500
    fs = exp(-s**2 * d)
    dfs = -2 * d * s * fs
    return fs, dfs


phi, dphi = f(y+0.2)
q0 = CF((0, -dphi))
mu0 = -dphi
wave.SetInitial(CF((q0, mu0)))
ng.Draw(u)
visoptions.scalfunction = 'u:3'
visoptions.subdivisions = '4'
visoptions.autoscale = False
visoptions.mmaxval = 8
visoptions.mminval = -8

t = 0
cnt = 0
dt = 0.1
redraw = 1

with ng.TaskManager():
    while t < 0.7:
        wave.Propagate()
        t += dt
        cnt += 1

        # input("continue?")

        if cnt % redraw == 0:
            print("{:.3f}".format(t))
            ng.Redraw(True)
