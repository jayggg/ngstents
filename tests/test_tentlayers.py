from ngsolve import Mesh
from ngstents import TentSlab
from netgen.geom2d import unit_square


def test_tent_nlayers():
    mesh = Mesh(unit_square.GenerateMesh(maxh=.3))
    dt = 5
    c = 1
    global_ct = 0.999
    method = "vol"
    heapsize = 5 * 1000 * 1000
    tentslab = TentSlab(mesh, method, heapsize)
    tentslab.SetMaxWavespeed(c)
    tentslab.PitchTents(dt, global_ct=global_ct)
    tents = [tentslab.GetTent(i) for i in range(tentslab.GetNTents())]
    layers = set([t.level for t in tents])
    assert len(layers) == tentslab.GetNLayers(), "Incorrect number of layers"
