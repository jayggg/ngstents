from netgen.geom2d import SplineGeometry
from ngsolve import (Mesh, L2, GridFunction, CoefficientFunction,
                     exp, x, y, TaskManager)
from ngstents import TentSlab
from ngstents.conslaw import Burgers


def test_burgers2D():
    order = 2
    dt = 0.025
    c = 16
    tend = dt
    method = "edge"

    geom = SplineGeometry()
    geom.AddRectangle((0, 0), (1, 1), bc=1)
    mesh = Mesh(geom.GenerateMesh(maxh=0.2))
    ts = TentSlab(mesh, method)
    ts.SetMaxWavespeed(c)
    success = ts.PitchTents(dt)
    try:
        assert success is True
    except AssertionError as e:
        msg = "Slab could not be pitched"
        e.args += ("Failed to pitch slab", msg)
        raise
    cf = CoefficientFunction(exp(-50*((x-0.3)*(x-0.3)+(y-0.3)*(y-0.3))))

    V = L2(mesh, order=order)
    u = GridFunction(V, "u")
    burg = Burgers(u, ts)
    burg.SetTentSolver("SARK", substeps=order*order)
    burg.SetInitial(cf)

    t = 0
    with TaskManager():
        while t < tend - dt/2:
            burg.Propagate()
            t += dt
