from netgen.geom2d import unit_square
from ngsolve import Mesh
from ngstents import TentSlab
from pytest import approx


def test_tent_properties():
    mesh = Mesh(unit_square.GenerateMesh(maxh=.2))
    dt = 0.05
    c = 16

    # Tent slab tests
    tentslab = TentSlab(mesh, dt, c)
    ntents = tentslab.GetNTents()
    assert ntents == 168
    slabheight = tentslab.GetSlabHeight()
    maxslope = tentslab.MaxSlope()
    assert maxslope < 0.1
    assert slabheight == 0.05

    tentslab.DrawPitchedTentsVTK()  # should create output.vtk

    # for ngsgui / tents_visualization and webgui
    results = tentslab.DrawPitchedTentsGL()
    tentdata, tenttimes, ntents, nlevels = results
    assert ntents == 168
    assert nlevels == 20
    assert len(tentdata) == 2872
    assert len(tenttimes) == 2872
    assert tentdata[0] == 0  # tent number
    assert tentdata[1] == 0  # level
    assert tentdata[2] == 0  # vertex number
    assert tentdata[3] == 0  # elnr
    assert tenttimes[3] == approx(0.0125)

    # Tent tests
    tent = tentslab.GetTent(0)
    assert tent.vertex == 0
    assert tent.ttop == approx(0.0125)
    assert tent.tbot == approx(0)
    nbv = tent.nbv.NumPy()
    nbtime = tent.nbtime.NumPy()
    els = tent.els.NumPy()
    internal_facets = tent.internal_facets.NumPy()
    assert len(nbv) == 2
    assert nbv[0] == 4          # vertex numbers
    assert nbv[1] == 19
    assert len(nbtime) == 2
    assert nbtime[0] == 0.0
    assert nbtime[1] == 0.0
    assert len(els) == 1
    assert els[0] == 0          # element number
    assert len(internal_facets) == 2
    assert internal_facets[0] == 0
    assert internal_facets[1] == 1
