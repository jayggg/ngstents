from netgen.meshing import Mesh, MeshPoint, Pnt, Element1D, Element0D
from ngsolve import Mesh as ngMesh
from ngstents import TentSlab


def test_tent_properties():
    nmesh = Mesh(dim=1)
    N = 10
    pnums = []
    for i in range(0, N+1):
        pnums.append(nmesh.Add(MeshPoint(Pnt(0.5**i, 0, 0))))
    idx = nmesh.AddRegion("mat", dim=1)
    for i in range(0, N):
        nmesh.Add(Element1D([pnums[i], pnums[i+1]], index=idx))
    idx_l = nmesh.AddRegion("left", dim=0)
    idx_r = nmesh.AddRegion("right", dim=0)
    nmesh.Add(Element0D(pnums[0], index=idx_l))
    nmesh.Add(Element0D(pnums[N], index=idx_r))
    mesh = ngMesh(nmesh)
    dt = 10
    c = 1
    method = "vol"
    # Tent slab tests
    tentslab = TentSlab(mesh, method,
                        5000000)
    tentslab.SetWavespeed(c)
    success = tentslab.PitchTents(dt)
    try:
        assert success is True
    except AssertionError as e:
        msg = "Slab could not be pitched"
        e.args += ("Failed to pitch slab", msg)
        raise
    ntents = tentslab.GetNTents()
    maxslope = tentslab.MaxSlope()
    try:
        assert maxslope <= 1.0/c
    except AssertionError as e:
        msg = "max slope = " + str(maxslope)
        msg += "   1.0/c = " + str(1.0/c)
        e.args += ("Failed at MaxSlopeTest!", msg)
        raise
    print(ntents)


if __name__ == "__main__":
    test_tent_properties()
