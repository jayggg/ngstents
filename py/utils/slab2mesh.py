from time import time
from collections import defaultdict

from netgen.geom2d import unit_square, SplineGeometry
import netgen.meshing as nm
import ngsolve as ng
import numpy as np

from ngstents import TentSlab
from ngstents.utils import Make1DMesh


class STv:
    """
    A simple class for a spacetime vertex

    Properties

    sv:    (ng.MeshNode) spatial vertex MeshNode
    front: (int) front number (0 is spatial mesh)
    pt:    (list) [x,y,z]
    stv:   (int) spacetime vertex nr (1-based)
    bnd:   (bool) True for a spacteime boundary vertex
    """

    def __init__(self, ngmesh, sv, front, scale=1.0, time=0.0,
                 dt=None, bvs=None):
        """
        ngmesh: (netgen.meshing.mesh)
        sv:     (int) spatial vertex NodeId
        front:  (int) front number (0 is spatial mesh)
        scale:  (float) scaling factor for time dimension
        time:   (float) central vertex time
        dt:     (float) final time for slab
        bvs:    (set) spatial boundary vertices
        """
        self.sv = sv
        self.front = front
        self.bnd = True if dt is None else sv.nr in bvs or time == dt
        self.pt = list(sv.point) + [time*scale]
        if len(self.pt) == 2:
            self.pt.append(0)
        self.stv = ngmesh.Add(nm.MeshPoint(nm.Pnt(*self.pt)))

    def __repr__(self):
        return "STv({},{},{},{},{})".format(
            self.sv, self.front, self.pt, self.stv, self.bnd)


class SlabConverter:
    """
    A class for converting a tent slab to a valid NGSolve mesh.
    Provides a method 'Convert' to generate the mesh.
    Once the conversion has been performed, a number of
    data structures are available as attributes of the converter object.
    These include 'mesh', which is the generated mesh and 'gfixmap',
    used with an ngstents conservation law to view time slices of a solution.
    'gfixmap' is a list of dicts, one for each front, mapping a tent vertex
    pitched in that front to the corresponding vertex of the generated mesh.
    If the constructor is called with p_hd = 2, the tent's internal edges
    are also mapped to edges of the generated mesh.
    """

    def __init__(self, tps, p_hd=1):
        """
        INPUTS

        tps: An ngstents N-D TentSlab instance, where N is one plus
        the dimension of the spatial mesh of the tent-pitched slab

        p_hd: The order of the high dimensional H1 space (1 or 2)
        """
        self.spatialmesh = tps.mesh
        self.N = tps.mesh.dim + 1
        self.tps = tps       # tent slab
        self.p_hd = p_hd     # order of the high dimensional H1 space
        if p_hd > 2:
            raise NotImplementedError("p_hd must be 1 or 2")
        self.dt = tps.GetSlabHeight()
        self.ntents = tps.GetNTents()
        self.nlayers = tps.GetNLayers()
        self.nfronts = self.nlayers + 1
        self.vertices = self.spatialmesh.vertices
        self.ngmesh = nm.Mesh(dim=self.N)   # ngmesh under construction
        self.tscale = 1.0    # time scaling factor
        self.vdata = {}      # dict {(front, vtx): STv}
        # list of spacetime points in vertex order, each pt a list of 3 floats
        self.stpts = []
        self.surfels = []    # list of surface elements (ccw vertex tuples)
        self.f2vs = []       # nested list, where outer index is front
        self.v2fs = []       # nested list, where outer index is vertex
        self.bverts = None   # set of spacetime boundary vertices
        self.mesh = None     # ngsolve N-D mesh
        self.gfixmap = None  # dict {(front, vtx): index into st vertices}

    def Convert(self, tscale=1.0):
        """
        Convert the tent pitched slab to an N-D mesh, providing
        timing information and counts.
        """
        self.tscale = tscale
        begin = start = time()
        self._AddVertices()
        print("add vertices {:.5f}".format(time()-start))
        start = time()
        self._AddVolumeElements()
        print("add volume elements {:.5f}".format(time()-start))
        start = time()
        self._AddSurfaceElements()
        print("add surface elements {:.5f}".format(time()-start))
        start = time()
        self.mesh = ng.Mesh(self.ngmesh)
        print("make ngsolve mesh {:.5f}".format(time()-start))
        start = time()
        self._MakeMap()
        print("make index map", time()-start)
        print("{} verts, {} vol elems, {} surf elems in {:.5f}.".format(
            self.mesh.nv, self.mesh.GetNE(ng.VOL), self.mesh.GetNE(ng.BND),
            time()-begin))

    def _tomeshv(self, vnr):
        """
        Convert a vertex number (int) to a MeshNode
        """
        return self.spatialmesh[ng.NodeId(ng.VERTEX, vnr)]

    def _AddVertices(self):
        """
        Add all 3D vertices to the mesh, first those associated with
        spatial vertices, then those associated with tent vertices.
        In the process, construct the dict vdata, the nested lists
        f2vs and v2fs, and the set bverts, used to construct elements.
        """
        mesh, spmesh = self.ngmesh, self.spatialmesh
        vdata, stpts = self.vdata, self.stpts
        bvs = {v.nr for el in spmesh.Elements(ng.BND)
               for v in spmesh[el].vertices}
        vs = []
        bverts = []
        for v in self.vertices:
            stv = STv(mesh, v, 0)
            vdata[(0, v.nr)] = stv
            stpts.append(stv.pt)
            self.v2fs.append([0])
            vs.append(v)
            bverts.append(stv.stv)
        self.f2vs = [vs] + [[] for i in range(self.nlayers)]

        for i in range(self.ntents):
            t = self.tps.GetTent(i)
            ft = t.level + 1
            v = self._tomeshv(t.vertex)
            stv = STv(mesh, v, ft, self.tscale, t.ttop, self.dt, bvs)
            vdata[(ft, t.vertex)] = stv
            stpts.append(stv.pt)
            self.v2fs[t.vertex].append(ft)
            self.f2vs[ft].append(v)
            if stv.bnd:
                bverts.append(stv.stv)
        self.bverts = set(bverts)

    def _AddVolumeElements(self):
        """
        For each front after the first, for each pitch vertex in the front,
        for each element incident on the pitch vertex, form an N-D element,
        adding it to the mesh.  We find each other vertex for such an
        N-D element by identifying its most recent prior front.
        """
        mesh = self.spatialmesh
        vdata, v2fs, f2vs = self.vdata, self.v2fs, self.f2vs
        ElementND = nm.Element2D if self.N == 2 else nm.Element3D
        idx_dom = self.ngmesh.AddRegion("mat", dim=self.N)
        for i, vs in enumerate(f2vs[1:]):
            ft = i+1
            def prevft(ev): return max(f for f in v2fs[ev] if f < ft)
            for v in vs:
                vdat = vdata[(ft, v.nr)]
                for el in v.elements:
                    elpnums = [vdat.stv] + \
                              [vdata[(prevft(elv.nr), elv.nr)].stv
                               for elv in mesh[el].vertices]
                    self.ngmesh.Add(ElementND(idx_dom, [*elpnums]))

    def _facetvs(self, vs):
        """
        Get the N+1 facets as lists of N vertices from the
        volume element's N+1 vertices
        """
        dim = self.N
        tvs = vs + vs[:dim-1]
        return [tvs[i:i+dim] for i in range(dim+1)]

    def _AddSurfaceElements(self):
        """
        Generate element facet tuples whose vertices are all surface vertices.
        Ensure correct orientation, then add them to the mesh.
        """
        bverts = self.bverts
        ElementsND = (self.ngmesh.Elements2D if self.N == 2
                      else self.ngmesh.Elements3D)
        idx_surf = self.ngmesh.AddRegion("surf", dim=self.N-1)
        for i, elNd in enumerate(ElementsND()):
            vlist = [v.nr for v in elNd.vertices]
            facets = self._facetvs(vlist)
            for facet in facets:
                # ensure that all vertices are on the surface
                if len(set(facet)-bverts) == 0:
                    self.surfels.append(facet)
        self._OrientSurfaceElements()
        for elt in self.surfels:
            if self.N == 2:
                el = nm.Element1D([*elt], index=idx_surf)
            else:
                el = nm.Element2D(idx_surf, elt)
            self.ngmesh.Add(el)

    def _OrientSurfaceElements(self):
        """
        Given a list of n surface elements, represented as tuples
        of integers, construct a corresponding n x N x N array of
        their spacetime point values and get the center point.
        Then orient the elements ccw by ensuring a positive
        determinant of a matrix with rows:
        - the vector from the center point to the first element point
        - the vectors formed by subtracting the first point from the
          other two points.
        These determinants are computed as vectorized scalar triple products.
        """
        # Note: pnums in surfel are 1-based
        pts = np.array([[self.stpts[v-1] for v in surfel]
                        for surfel in self.surfels])
        p = pts[:, 0, :] - np.mean(pts, axis=(0, 1))
        if self.N == 2:
            sl = pts[:, 0, :][:, np.newaxis, :]
            v = (pts - sl)[:, 1, :]
            test = np.cross(p, v)[:, 2]
        else:
            sl = pts[:, 0, :][:, np.newaxis, :]
            pts = (pts - sl)
            v, w = pts[:, 1, :], pts[:, 2, :]
            test = np.einsum('ij, ij->i', p, np.cross(v, w))
        test = test < 0
        self.surfels = [elt if test[i] else [elt[1], elt[0]] + elt[2:]
                        for i, elt in enumerate(self.surfels)]

    def _MakeMap(self):
        smesh = self.spatialmesh
        mesh = self.mesh
        v2fs = self.v2fs
        f2vs = self.f2vs
        nv = mesh.nv
        nvs = smesh.nv

        tmap = defaultdict(dict)
        for k, v in self.vdata.items():
            tmap[k[0]][k[1]] = v.stv.nr - 1
            gfixmap = [tmap[f] for f in range(self.nfronts)]
        if self.p_hd > 1:
            # *Temporarily* work around issue:
            # https://ngsolve.org/forum/ngspy-forum/
            #   1413-original-edges-retained-after-mesh-refine#3812
            Z = ng.HCurl(smesh, order=0)
            fd = Z.FreeDofs()
            true_edges = [e for i, e in enumerate(smesh.edges) if fd[i]]
             
            minedgenr = min(e.nr for e in true_edges)
            # extend gfixmap to include edge dof indices
            # Note that edge.vertices always returns a pair
            # but vertex.edges includes edges not incident on the vertex.

            # first construct a dict from ordered pair of st vertices to st edge
            e2vs = {e.nr: [v.nr for v in mesh[e].vertices]
                    for e in mesh.edges}
            vs2e = {tuple(v): k for k, v in e2vs.items()}
            #print('vs2e', vs2e)

            # for front 0, map each spatial edge to the corresponding
            # edge of the spacetime mesh, then use the offsets nvs-minedgenr
            # and nv to map spatial dof nrs for front 0 to spacetime dof nrs.
            edict = {}
            for e in true_edges:
                svs = [v.nr for v in smesh[e].vertices]
                stvs = [gfixmap[0][v.nr] for v in smesh[e].vertices]
                # print("stvs", stvs)
                gfixmap[0][nvs + e.nr - minedgenr] = nv + vs2e[tuple(stvs)]

            # for each front > 0, for each spatial vertex in the front
            # get the corresponding spacetime vertex and
            # the spatial vertex' edges.
            # for each of these edges, map the opposite vertex
            # to the spacetime vertex and look up the spacetime edge,
            # then use the offsets nvs-minedgenr and nv to map spatial dof
            # nrs for the front to spacetime dof nrs.
            # Mapping the vertices to spacetime vertices requires that we find
            # the nearest previous front
            for i, vs in enumerate(f2vs[1:]):
                ft = i+1
                for v in self.f2vs[ft]:
                    stv = gfixmap[ft][v.nr]
                    for e in smesh[v].edges:
                        sev = [ev for ev in smesh[e].vertices if ev != v][0]
                        evf = max(f for f in v2fs[sev.nr] if f < ft)
                        stevs = tuple(sorted([stv, gfixmap[evf][sev.nr]]))
                        ste = vs2e[stevs]
                        gfixmap[ft][nvs + e.nr - minedgenr] = nv + ste
        self.gfixmap = gfixmap


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convenience functions for testing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def Diskmesh(maxh=.1):
    geo = SplineGeometry()
    geo.AddCircle([0, 0], 1)
    return geo.GenerateMesh(maxh=maxh)


def Ringmesh(maxh=.1):
    geo = SplineGeometry()
    geo.AddCircle([0, 0], 1, leftdomain=0, rightdomain=1)
    geo.AddCircle([0, 0], 2)
    return geo.GenerateMesh(maxh=maxh)


def Trigmesh(maxh=.1):
    geo = SplineGeometry()
    vs = [[0, 0], [1, 0], [.5, .7071]]
    pts = [geo.AddPoint(*v) for v in vs]
    _ = [geo.Append([pts[i], pts[(i+1) % 3]]) for i in [0, 1, 2]]
    return geo.GenerateMesh(maxh=maxh)


if __name__ == '__main__':
    # SIMULATION SETTINGS
    dim = 3
    dt = 0.09
    local_ctau = True
    global_ctau = .999
    wavespeed = 6
    tscale = 5.0
    # settings used only when dim = 3
    geotype = "ring"
    maxh = .2
    order_hd = 2

    import netgen.gui
    if dim == 2:
        mesh = ng.Mesh(Make1DMesh([[0, 1]], [10]))
    else:
        if geotype == "square":
            msh = unit_square.GenerateMesh(maxh=maxh)
        elif geotype == "trig":
            msh = Trigmesh(maxh=maxh)
        elif geotype == "disk":
            msh = Diskmesh(maxh=maxh)
        elif geotype == "ring":
            msh = Ringmesh(maxh=maxh)
        else:
            raise ValueError("Undefined geotype")
        mesh = ng.Mesh(msh)
        mesh.Refine()

    ts = TentSlab(mesh, method="edge")
    ts.SetMaxWavespeed(wavespeed)
    ts.PitchTents(dt=dt, local_ct=local_ctau, global_ct=global_ctau)
    sc = SlabConverter(ts, order_hd)
    sc.Convert(tscale=tscale)

    fes = ng.H1(sc.mesh, order=order_hd)
    gf = ng.GridFunction(fes)
    gf.vec[:] = 0
    ix = sc.gfixmap[0][0]
    gf.vec.FV()[ix] = 1.0
    nv = len(sc.vertices)
    ix = sc.gfixmap[0][nv]
    gf.vec.FV()[ix] = 5.0
    ix = sc.gfixmap[0][nv+1]
    gf.vec.FV()[ix] = 5.0
    ng.Draw(gf)
