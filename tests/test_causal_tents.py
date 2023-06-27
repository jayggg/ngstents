"""
Module test_causal_tents

Tests edge-based and volume-based methods in 1D, 2D and 3D

Note: In all tests, setting the global_ct slightly less than 1 ensures
      a causal slab.  If global_ct=1, the default, then roundoff error may
      cause tests to fail.
"""
from ngsolve import Mesh
from ngstents import TentSlab
from ngsolve.meshes import Make1DMesh


def Get1DMesh():
    return Make1DMesh(10)


def Get2DMesh():
    from netgen.geom2d import unit_square
    mesh = Mesh(unit_square.GenerateMesh(maxh=.3))
    mesh.Refine()
    return mesh


def Get3DMesh():
    from netgen.csg import unit_cube
    mesh = Mesh(unit_cube.GenerateMesh(maxh=.5))
    mesh.Refine()
    return mesh


def test_1D_vol_causal():
    mesh = Get1DMesh()
    dt = 10
    c = 1
    global_ct = 0.999
    method = "vol"
    heapsize = 5 * 1000 * 1000
    tentslab = TentSlab(mesh, method, heapsize)
    tentslab.SetMaxWavespeed(c)
    success = tentslab.PitchTents(dt, global_ct=global_ct)
    assert success, "Slab could not be pitched"
    maxslope = tentslab.MaxSlope()
    expected = 1.0 / c
    msg = "max slope {} exceeded {}".format(maxslope, expected)
    assert maxslope <= expected, msg


def test_1D_edge_causal():
    mesh = Get1DMesh()
    dt = 10
    c = 1
    global_ct = 0.999
    method = "edge"
    heapsize = 5 * 1000 * 1000
    tentslab = TentSlab(mesh, method, heapsize)
    tentslab.SetMaxWavespeed(c)
    success = tentslab.PitchTents(dt, global_ct=global_ct)
    assert success, "Slab could not be pitched"
    maxslope = tentslab.MaxSlope()
    expected = 1.0 / c
    msg = "max slope {} exceeded {}".format(maxslope, expected)
    assert maxslope <= expected, msg


def test_2D_vol_causal():
    mesh = Get2DMesh()
    dt = 10
    c = 1
    global_ct = 0.999
    method = "vol"
    heapsize = 5 * 1000 * 1000
    tentslab = TentSlab(mesh, method, heapsize)
    tentslab.SetMaxWavespeed(c)
    success = tentslab.PitchTents(dt, local_ct=True, global_ct=global_ct)
    assert success, "Slab could not be pitched"
    maxslope = tentslab.MaxSlope()
    expected = 1.0 / c
    msg = "max slope {} exceeded {}".format(maxslope, expected)
    assert maxslope <= expected, msg


def test_2D_edge_causal():
    mesh = Get2DMesh()
    dt = 10
    c = 1
    global_ct = 0.999
    method = "edge"
    heapsize = 5 * 1000 * 1000
    tentslab = TentSlab(mesh, method, heapsize)
    tentslab.SetMaxWavespeed(c)
    success = tentslab.PitchTents(dt, local_ct=True, global_ct=global_ct)
    assert success, "Slab could not be pitched"
    maxslope = tentslab.MaxSlope()
    expected = 1.0 / c
    msg = "max slope {} exceeded {}".format(maxslope, expected)
    assert maxslope <= expected, msg


def test_3D_vol_causal():
    mesh = Get2DMesh()
    dt = 10
    c = 1
    global_ct = 0.999
    method = "vol"
    heapsize = 5 * 1000 * 1000
    tentslab = TentSlab(mesh, method, heapsize)
    tentslab.SetMaxWavespeed(c)
    success = tentslab.PitchTents(dt, local_ct=True, global_ct=global_ct)
    assert success, "Slab could not be pitched"
    maxslope = tentslab.MaxSlope()
    expected = 1.0 / c
    msg = "max slope {} exceeded {}".format(maxslope, expected)
    assert maxslope <= expected, msg


def test_3D_edge_causal():
    mesh = Get2DMesh()
    dt = 10
    c = 1
    global_ct = 0.999
    method = "edge"
    heapsize = 5 * 1000 * 1000
    tentslab = TentSlab(mesh, method, heapsize)
    tentslab.SetMaxWavespeed(c)
    success = tentslab.PitchTents(dt, local_ct=True, global_ct=global_ct)
    assert success, "Slab could not be pitched"
    maxslope = tentslab.MaxSlope()
    expected = 1.0 / c
    msg = "max slope {} exceeded {}".format(maxslope, expected)
    assert maxslope <= expected, msg
