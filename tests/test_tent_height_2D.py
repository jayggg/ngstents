from netgen.geom2d import unit_square
from ngsolve import Mesh, NodeId, FACET, VERTEX
from ngstents import TentSlab
import numpy as np


def test_tent_height():
    """Partial test for tents' height

    Check if the tents' height are not bigger than any of the neighbouring
    edges. Passing this test does NOT imply the fulfillment of causality
    conditions."""

    mesh = Mesh(unit_square.GenerateMesh(maxh=.2))
    dt = 0.05
    c = 16
    tol = 1e-12
    # Tent slab tests
    tentslab = TentSlab(mesh, dt, c)
    ntents = tentslab.GetNTents()
    for itent in range(ntents):
        tent = tentslab.GetTent(itent)
        tent_v = tent.vertex
        tent_pt = mesh[NodeId(VERTEX, tent_v)].point
        time_center = tent.ttop
        tent_facets = tent.internal_facets
        for edge_nr in tent_facets:
            edge = mesh[NodeId(FACET, edge_nr)]
            edge_v = edge.vertices[1].nr \
                if edge.vertices[0].nr == tent_v  \
                else edge.vertices[0].nr
            other_v = edge.vertices[0].nr \
                if edge.vertices[0].nr == tent_v \
                else edge.vertices[1].nr
            try:
                assert (tent_v == other_v)
            except AssertionError as e:
                msg = "edge_v = "+str(edge_v)+" tent_v = "+str(tent_v)
                msg += " edge.v[0] = " + str(edge.vertices[0].nr)
                msg += " edge.v[1] = " + str(edge.vertices[1].nr)
                e.args += ("ERROR on tent data structure", msg)
                raise
            edge_pt = mesh[NodeId(VERTEX, edge_v)].point
            edge_pt_local_id = np.where(tent.nbv.NumPy() == edge_v)[0][0]
            dist = np.linalg.norm(np.subtract(tent_pt, edge_pt))
            time_node = tent.nbtime[edge_pt_local_id]
            diff_t = time_center - time_node

            try:

                assert (dist/c >= diff_t - tol)

            except AssertionError as e:
                msg = "tent id = " + str(itent)+" tent_v = " + str(tent_v)
                msg += " edge_v = " + str(edge_v) + \
                    " time_center = " + str(time_center)
                msg += " time_node = "+str(tent.nbtime[edge_pt_local_id])
                msg += " coord center = " + str(tent_pt)
                msg += " coord node = " + str(edge_pt)
                msg += " edge length = " + str(dist)
                for iv in range(len(tent.nbtime)):
                    msg += " vertex = " + \
                        str(tent.nbv[iv]) + " time = " + str(tent.nbtime[iv])
                e.args += ("ERROR: tent has slope bigger than velocity" +
                           " along edge!", msg)
                raise
        try:

            assert (1.0/c >= tent.MaxSlope())

        except AssertionError as e:
            msg = "tent id = " + str(itent)+" tent_v = " + str(tent_v)
            msg += "max slope = " + str(tent.MaxSlope())
            msg += "1/c = " + str(1.0/c)
            e.args += ("ERROR: tent has slope bigger than velocity!", msg)
            raise


if __name__ == "__main__":
    test_tent_height()
