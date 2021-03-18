from netgen.geom2d import unit_square
from ngsolve import Mesh
from ngstents import TentSlab


def test_tent_properties():
    mesh = Mesh(unit_square.GenerateMesh(maxh=.2))
    mesh.Refine()
    mesh.Refine()

    dt = 0.5
    c = 1
    method = "edge"
    local_ctau = True
    # Tent slab tests
    tentslab = TentSlab(mesh, method, 1000000000)
    tentslab.SetMaxWavespeed(c)

    maxslope = 1.0/c
    global_ctau = 1
    while maxslope >= 1.0/c:
        print("Pitching with global cte = {}".format(global_ctau))
        success = tentslab.PitchTents(dt, local_ctau, global_ctau)
        try:
            assert success is True
        except AssertionError as e:
            msg = "Slab could not be pitched"
            e.args += ("Failed to pitch slab", msg)
        ntents = tentslab.GetNTents()
        maxslope = tentslab.MaxSlope()
        print("Pitched {} tents with max slope = {}".format(ntents, maxslope))
        global_ctau /= maxslope

    try:
        assert maxslope <= 1.0/c
    except AssertionError as e:
        msg = "max slope = " + str(maxslope)
        msg += "   1.0/c = " + str(1.0/c)
        e.args += ("Failed at MaxSlopeTest!", msg)


    # for ngsgui / tents_visualization and webgui
    results = tentslab.DrawPitchedTentsGL()
    tentdata, tenttimes, ntents, nlevels = results



if __name__ == "__main__":
    test_tent_properties()
