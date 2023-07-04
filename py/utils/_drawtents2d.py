from netgen.meshing import Mesh, FaceDescriptor, Element3D, Element2D, \
    MeshPoint, Pnt
import ngsolve as ng
import numpy as np


def DrawPitchedTents(self, uptolevel=None):  # self = TentSlab py object
    """
    Make a 2+1 dimensional mesh of spacetime tents, drawable using netgen.
    If uptolevel=L, then only tents of level L and lower are plotted,
    and if uptolevel is not given, all tents are included in the mesh.

    RETURNS:
      mesh3d: the spacetime mesh
      vertexlevels: tentlevels assigned to apex vertices (H1 function)
      tetlevels: tentlevels assigned to spacetime tetrahedra (DG function)
    """

    if self.mesh.dim != 2:
        raise NotImplementedError("Only supported for 2D spatial mesh")
    else:

        if self.mesh.ngmesh.GetNrIdentifications():
            raise NotImplementedError(
                "Don't know a good way to show periodic tents")

        # Make new a 3D mesh
        m = Mesh()
        fdbot = m.Add(FaceDescriptor(bc=1, domin=1, surfnr=1))
        fdwal = m.Add(FaceDescriptor(bc=2, domin=1, surfnr=2))
        fdtop = m.Add(FaceDescriptor(bc=3, domin=1, surfnr=3))
        bottomspatialmesh = self.mesh.ngmesh
        pmap = {}
        levelsnodal = {}

        for e in bottomspatialmesh.Elements2D():
            for v in e.vertices:
                if (v not in pmap):
                    pmap[v] = m.Add(bottomspatialmesh[v])
                    levelsnodal[v] = 0

        for e in bottomspatialmesh.Elements2D():
            m.Add(Element2D(fdbot, [pmap[v] for v in e.vertices]))

        # make:  tentsbylayer[layer] = list of tent nums in layer
        tentsbylayers = {i: [] for i in range(self.GetNLayers())}
        for i in range(self.GetNTents()):
            t = self.GetTent(i)
            tentsbylayers[t.level].append(i)

        bottomspatialelements = list(bottomspatialmesh.Elements2D())
        front = pmap
        levelstet = {}

        if uptolevel is None:
            uptolevel = self.GetNLayers()

        # Add surface elements and volume elements layer by layer
        for layer in range(uptolevel):
            for tn in tentsbylayers[layer]:
                t = self.GetTent(tn)
                pxy = self.mesh.vertices[t.vertex].point
                vtop = m.Add(MeshPoint(Pnt(pxy[0], pxy[1], t.ttop)))
                levelsnodal[vtop] = layer + 1
                vbot = front[t.vertex + 1]
                front[t.vertex + 1] = vtop

                # Make vertical (const in space) boundary elements
                for e in self.mesh.Elements(ng.BND):
                    if e.vertices[0].nr == t.vertex:
                        o = front[e.vertices[1].nr + 1]
                        m.Add(Element2D(fdwal, [vtop, vbot, o]))
                        # implicitly used orientation of e in above ordering
                    if e.vertices[1].nr == t.vertex:
                        o = front[e.vertices[0].nr + 1]
                        m.Add(Element2D(fdwal, [vtop, o, vbot]))

                # Make advancing front faces and new elements underneath
                for el in t.els:
                    topface = [
                        front[v] for v in bottomspatialelements[el].vertices
                    ]
                    m.Add(Element2D(fdtop, topface))
                    tet = m.Add(Element3D(1, topface + [vbot]))
                    levelstet[tet] = layer + 1

        m.Update()
        mesh3d = ng.Mesh(m)
        W = ng.L2(mesh3d, order=0)
        tetlevels = ng.GridFunction(W)
        tetlevels.vec.FV().NumPy()[:] = np.array(
            [levelstet[tet] for tet in levelstet])
        V = ng.H1(mesh3d, order=1)
        vertexlevels = ng.GridFunction(V)
        vertexlevels.vec.FV().NumPy()[:] = np.array(
            [levelsnodal[v] for v in levelsnodal])

        return mesh3d, vertexlevels, tetlevels
