from netgen.geom2d import SplineGeometry
from ngsolve import Mesh, CoefficientFunction, exp, x, y, TaskManager
from ngstents import Burgers, TentSlab


def test_burgers2D():
    order = 2
    dt = 0.025
    c = 16
    tend = dt

    geom = SplineGeometry()
    geom.AddRectangle((0, 0), (1, 1), bc=1)
    mesh = Mesh(geom.GenerateMesh(maxh=0.2))
    ts = TentSlab(mesh, dt, c)
    cf = CoefficientFunction(exp(-50*((x-0.3)*(x-0.3)+(y-0.3)*(y-0.3))))

    burg = Burgers(ts, order=order)
    sol = burg.sol
    burg.SetInitial(cf)

    t = 0
    with TaskManager():
        while t < tend - dt/2:
            burg.PropagatePicard(sol.vec, steps=order*order)
            t += dt