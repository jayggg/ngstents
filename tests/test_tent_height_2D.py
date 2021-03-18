from netgen.geom2d import unit_square
from ngsolve import Mesh, NodeId, FACET, VERTEX
from ngstents import TentSlab
import numpy as np


def test_tent_height():
    """
    Check that the edge slopes do not exceed 1/wavespeed.
    Passing this test does NOT imply the fulfillment of causality
    conditions.
    """

    mesh = Mesh(unit_square.GenerateMesh(maxh=.5))
    nref = 2
    for i in range(nref):
        mesh.Refine()
    dt = 0.05
    c = 16
    tol = 1e-12
    method = "vol"
    # Tent slab tests
    tentslab = TentSlab(mesh, method, 10**7)
    tentslab.SetMaxWavespeed(c)
    assert tentslab.PitchTents(dt), "Slab could not be pitched"
    for i in range(tentslab.GetNTents()):
        tent = tentslab.GetTent(i)
        tent_v = tent.vertex
        tent_pt = mesh[NodeId(VERTEX, tent_v)].point
        time_center = tent.ttop
        tent_facets = tent.internal_facets
        for edge_nr in tent_facets:
            edge = mesh[NodeId(FACET, edge_nr)]
            if edge.vertices[1].nr == tent_v:
                edge_v, other_v = edge.vertices[0].nr, edge.vertices[1].nr
            else:
                edge_v, other_v = edge.vertices[1].nr, edge.vertices[0].nr
            assert tent_v == other_v, "Error in tent structure"
            edge_pt = mesh[NodeId(VERTEX, edge_v)].point
            edge_pt_local_id = np.where(np.array(tent.nbv) == edge_v)[0][0]
            dist = np.linalg.norm(np.subtract(tent_pt, edge_pt))
            time_node = tent.nbtime[edge_pt_local_id]
            diff_t = time_center - time_node
            assert dist/c >= diff_t - tol, "Tent slope exceeds 1/c along edge!"
