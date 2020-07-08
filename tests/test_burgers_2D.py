from netgen.geom2d import SplineGeometry
from ngsolve import Mesh, CoefficientFunction, exp, x, y, TaskManager
from tents import ConsLaw


def test_conslaw_properties():
    order = 2
    dt = 0.025
    c = 16
    tend = dt

    geom = SplineGeometry()
    geom.AddRectangle((0, 0), (1, 1), bc=1)
    mesh = Mesh(geom.GenerateMesh(maxh=0.2))
    cf = CoefficientFunction(exp(-50*((x-0.3)*(x-0.3)+(y-0.3)*(y-0.3))))

    cl = ConsLaw(mesh, "burgers", order=order)
    sol = cl.sol
    cl.SetInitial(cf)
    cl.PitchTents(dt, c)
    cl.DrawPitchedTentsVTK('conslaw_tents')
    results = cl.DrawPitchedTentsGL()
    tentdata, tenttimes, ntents, nlevels = results

    assert cl.MaxSlope() < 0.1
    t = 0
    with TaskManager():
        while t < tend - dt/2:
            cl.PropagatePicard(sol.vec, steps=order*order)
            t += dt
