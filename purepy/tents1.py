"""
Module tents1

A pure Python implementation of a tent pitched slab with a number of
features and options for experimentation and debugging

Visualizing the tent-pitched slab:
For 1D or 3D trials, running in iPython allows visualization via pyplot
For 2D trials, running from the command line generates an HTML file.
"""

from collections import defaultdict

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import ngsolve as ng

from ngstents import Make1DMesh, Make1DMeshSpecified
from tentswebgui import Draw


class Tent(object):
    """
    A simple structure representng a space-time tent
    """

    def __init__(self, vertex, tbot, ttop, level, els):
        self.vertex = vertex
        self.tbot = tbot
        self.ttop = ttop
        self.level = level
        self.nbv = []
        self.nbtime = []
        self.els = els

    def __str__(self):
        nbv = ", ".join(["{}".format(nb) for nb in self.nbv])
        nbtime = ", ".join(["{:.3f}".format(nbt) for nbt in self.nbtime])
        s = "v: {}, tbot: {:.3f}, ttop: {:.3f}," + \
            " level: {},\n  nbts: {}\n  nbvs: {}"
        return s.format(self.vertex, self.tbot, self.ttop,
                        self.level, nbtime, nbv)


class TentPitchedSlab(object):
    """
    A collection of Tent objects whose union is a spacetime cylinder
    with base the spatial mesh and height the slab height provided
    to the pitch_tents method.
    """

    def __init__(self, mesh):
        """
        INPUTS

        mesh: an NGSolve mesh of any dimension
        """
        self.tents = []
        self.mesh = mesh
        self.dim = mesh.dim
        self.els = list(mesh.Elements(ng.VOL))
        self.edges = list(mesh.edges)
        self.vertices = list(mesh.vertices)
        self.v2e = {v: list(v.edges) for v in self.vertices}
        self.v2v = {v: list(set([v1 for e in self.v2e[v]
                                 for v1 in mesh[e].vertices])-set([v]))
                    for v in self.vertices}
        self.vpt = {v: mesh[v].point for v in self.vertices}
        self.slabheight = None
        self.wavespeed = None
        self.acceleration = None
        self.cbar = None
        self.tau = None
        self.edgedts = None
        self.tol = 1.0e-15
        self.printcheck = False
        self.printdisc = False
        self.maxerr = 0.
        self.maxgradtau = 0.
        self.edgelengths = None
        self.ve2qf = None
        self.ready_vs = None
        self.ready_factor = None
        self.reset_rf = False
        self.dt_method = None
        self.make_tent_method = None
        self.maxready = False
        self.qf = False
        self.edge_based = False
        self.Ctau = 1.0

    def print_tents(self):
        for i, t in enumerate(self.tents):
            print("tent ", i, t, "\n")

    def GetTent(self, tentnr):
        return self.tents[tentnr]

    def set_v2elgradphi(self):
        """
        construct a dict with key a vertex, value a dict with key
        elnr and value a gradient vector as a numpy array
        We can evaluate the gradient using a temp H1 space by getting
        an interior point and either calling CalcDShape on the element
        for the point or by setting each shape function using a
        GridFunction for the space and evaluating its grad.
        """
        fes = ng.H1(self.mesh, order=1)
        self.v2elgradphi = {}
        for vx in self.mesh.vertices:
            gradvals = {}
            for el in vx.elements:
                vs = list(self.mesh[el].vertices)
                pts = np.array([list(self.mesh[v].point) for v in vs])
                intpt = pts.mean(axis=0)
                fel = fes.GetFE(el)
                dshp = fel.CalcDShape(self.mesh(*intpt)).NumPy()
                vix = vs.index(vx)
                dshp = dshp[vix]
                if hasattr(dshp, "__iter__"):
                    gradvals[el.nr] = (np.array(dshp))
                else:  # 1D is special
                    gradvals[el.nr] = np.array([dshp])
            self.v2elgradphi[vx.nr] = gradvals

    def set_edgelengths(self):
        """
        Construct a dict that maps each tuple of adjacent vertices to
        the length of their connecting edge
        """
        def calcdist(vs):
            pts = np.array([self.vpt[v] for v in vs])
            return (((pts[0]-pts[1])**2).sum())**0.5

        self.edgelengths = {}
        for edge in self.edges:
            vs = list(edge.vertices)
            v0, v1 = vs
            length = calcdist(vs)
            self.edgelengths[(v0.nr, v1.nr)] = length
            self.edgelengths[(v1.nr, v0.nr)] = length

    def set_ve2qf(self):
        """
        construct a dict that maps a central vertex and an incident element
        to a QF_2 value.
        """
        mesh = self.mesh

        def getqf(vx, edges):
            adj = []
            for edge in edges:
                vs = list(mesh[edge].vertices)
                v0, v1 = vs
                if vx in vs:
                    adj.append(self.edgelengths[(v0.nr, v1.nr)])
                else:
                    opp = self.edgelengths[(v0.nr, v1.nr)]
            return min(1., opp/max(adj))

        self.set_edgelengths()
        self.ve2qf = {}
        for vx in self.mesh.vertices:
            for el in vx.elements:
                if self.dim == 2:
                    self.ve2qf[(vx.nr, el.nr)] = getqf(vx, mesh[el].edges)
                elif self.dim == 3:
                    qfs = []
                    for f in mesh[el].faces:
                        mf = mesh[f]
                        if vx in list(mf.vertices):
                            qfs.append(getqf(vx, mf.edges))
                    self.ve2qf[(vx.nr, el.nr)] = min(qfs)

    def pitch_tents(self, slabheight, wavespeed,
                    acceleration=0.0, checkgrad=False,
                    tol=1.0e-15, printcheck=False, printdisc=False,
                    ready_factor=0.5, reset_rf=False, maxready=False,
                    qf=False, edge_based=False, Ctau=1.0):
        """
        Generates all tents to fill the tent-pitched slab

        INPUTS

        slabheight: float - the time dimension of the slab to be generated
        wavespeed:  NGSolve CoefficientFunction - the initial wavespeed at
                      each spatial coordinate
        acceleration: optional float - a change in wavespeed to be applied
                      at each tent pitching step.
        checkgrad: optional boolean - print the computed grad norm, 1/cbar
                   and the error for each facet based on the dt value just
                   computed.
        tol: float (default 1.0e-14) - tolerance for gradient check.
        printcheck: optional boolean - if set, prints a line in gradient check
        ready_factor: initial ready factor
        maxready: optional boolean - if set, choose the ready vertex with
                  minimal layer with maximal ktilde.
        qf: optional boolean - if set, apply mesh quality factor to
            computed kbar values
        edge_based: if set, use edge_based algorithm
        Ctau: optional mesh constant for scaling ktilde used by the
               edge_based algorithm
        """
        self.slabheight = slabheight
        self.wavespeed = wavespeed
        self.acceleration = acceleration
        self.checkgrad = checkgrad
        self.tol = tol
        self.printcheck = printcheck
        self.printdisc = printdisc
        self.maxready = maxready
        self.ready_factor = ready_factor
        self.reset_rf = reset_rf
        self.qf = qf
        if self.qf:
            self.set_ve2qf()
        self.edge_based = edge_based
        self.dt_method = self.compute_dt_edgebased if edge_based \
            else self.compute_dt_H1
        self.Ctau = Ctau
        self.set_v2elgradphi()

        # initial advancing front at vertices
        tau = self.tau = {v: 0 for v in self.vertices}
        # wavespeed for each vertex
        self.cbar = {v: wavespeed(self.mesh(*self.vpt[v]))
                     for v in self.vertices}
        if min(self.cbar.values()) <= 1e-16:
            raise ValueError("Wavespeed is not strictly positive")

        if self.edge_based:
            self.edgedts = {e: self.compute_edge_dt(e)
                            for e in self.edges}
            vertdts = {v: max([self.edgedts[e] for e in v.edges])
                       for v in self.vertices}
        else:
            # reference heights for vertices for gauging progress potential
            vertdts = {v: self.dt_method(v) for v in self.vertices}

        # optimal heights for vertices at each iteration
        ktilde = vertdts.copy()

        # vertices for which good progress can be made
        self.ready_vs = self.vertices[:]

        # initial level for each vertex
        vs_level = {v: 0 for v in self.vertices}

        mintau, maxktilde = 0, max(ktilde.values())

        while mintau < slabheight and maxktilde > 0:
            # ensure that ready_vs is not an empty collection or return
            orig_rf = self.ready_factor
            while len(self.ready_vs) == 0 \
                    and maxktilde > 0:
                print("No ready vertices. Reducing ready_factor.")
                self.ready_factor /= 2
                for v in self.vertices:
                    if tau[v] < slabheight:
                        ktilde[v] = self.dt_method(v)
                        if ktilde[v] >= self.ready_factor * vertdts[v]:
                            self.ready_vs.append(v)
                maxktilde = max(ktilde.values())
            if self.reset_rf:
                self.ready_factor = orig_rf

            if len(self.ready_vs) == 0:
                return tau, ktilde

            # select a tent vertex and remove it from ready_vs
            minlevel = min(vs_level[v] for v in self.ready_vs)
            levelvs = [v for v in self.ready_vs if vs_level[v] == minlevel]
            if self.maxready:
                vi = sorted(levelvs, key=lambda v: ktilde[v])[-1]
            else:
                vi = levelvs[0]
            self.ready_vs.remove(vi)

            # pitch tent, update tau, ready_vs, ktilde if greedy method
            t = self.make_tent(
                vi, slabheight, vs_level, ktilde, vertdts)

            if checkgrad:
                self.checktent(t, tau, tol=self.tol)

            # update speeds (just for testing...)
            if self.acceleration != 0:
                self.cbar = {k: v + self.acceleration
                             for k, v in self.cbar.items()}
                # since we updated speeds, we need to update all the ktildes
                ktilde = {v: self.dt_method(v) for v in self.vertices}
                self.ktilde = ktilde
                self.ready_vs = [v for v in self.vertices
                                 if ktilde[v] > self.ready_factor*vertdts[v]
                                 and tau[v] < slabheight]

            mintau = min(tau.values())
            maxktilde = max(ktilde.values())
        return tau, ktilde

    def make_tent(self, vi, slabheight, vs_level, ktilde, vertdts):
        """
        Generate a tent using the greedy algorithm
        Update ktilde for neighbor vertices, then depending on the value,
        possibly update ready_vs.
        """
        tau = self.tau
        # construct the tent and append it
        t = Tent(vi, tau[vi],
                 min(slabheight, tau[vi]+ktilde[vi]), vs_level[vi],
                 [e.nr for e in self.mesh[vi].elements])
        self.tents.append(t)
        # update the advancing front
        tau[vi] = t.ttop

        # update the level for the central vertex
        vs_level[vi] += 1

        # keep ktilde consistent
        ktilde[vi] = 0

        for nb in self.v2v[vi]:
            t.nbv.append(nb.nr)
            t.nbtime.append(tau[nb])
            # update the levels for neighbor vertices
            if vs_level[nb] < t.level + 1:
                vs_level[nb] = t.level + 1

            # update the ktilde values for neighbor vertices
            if self.tau[nb] >= self.slabheight:
                continue
            dtval = self.dt_method(nb)
            ktilde[nb] = dtval
            if dtval > self.ready_factor*vertdts[nb] \
                    and nb not in self.ready_vs:
                self.ready_vs.append(nb)
            # prevent degenerate tents on refined meshes
            elif dtval < self.ready_factor*vertdts[nb] \
                    and nb in self.ready_vs:
                self.ready_vs.remove(nb)

        return t

    def compute_dt_edgebased(self, vtx, ret_dts=False):
        """
        Compute the dt for a vertex
        """
        tau = self.tau
        dts = [tau[v] - tau[vtx] + self.edgedts[e]
               for v in self.v2v[vtx] for e in self.v2e[v]]
        if ret_dts:
            return dts
        else:
            return max(0, min(dts))

    def compute_edge_dt(self, edge):
        """
        Compute the reference dt for an edge from its length and speed
        """
        speed = max(self.cbar[v] for v in edge.vertices)
        pts = np.array([self.vpt[v] for v in edge.vertices])
        diffs = pts[0, :] - pts[1, :]
        edge_length = ((diffs**2).sum())**0.5
        return edge_length / speed * self.Ctau

    def compute_dt_H1(self, vtx, ret_dts=False):
        """
        Compute the dt for a vertex
        Optionally return instead the list of dt values
        computed for each facet opposite to the vertex.
        """
        dim = self.dim
        dts = []
        for el in self.mesh[vtx].elements:
            verts = list(self.mesh[el].vertices)
            verts.remove(vtx)

            cbar = max(self.cbar[v] for v in verts + [vtx])
            t0 = self.tau[vtx]
            gphi0 = self.v2elgradphi[vtx.nr][el.nr]
            a = sum(gphi0**2)
            b = 2 * sum(self.tau[v] * self.v2elgradphi[v.nr][el.nr].dot(gphi0)
                        for v in verts)
            c = sum(self.tau[v] * self.tau[v2] *
                    self.v2elgradphi[v.nr][el.nr].dot(
                        self.v2elgradphi[v2.nr][el.nr])
                    for v in verts for v2 in verts) - 1./cbar**2
            # normalize to minimize roundoff errors
            a, b, c = 1., b/a, c/a
            disc = b**2 - 4*a*c
            if disc >= 0:
                newt0 = (np.sqrt(disc) - b)/2/a
            else:
                if self.printdisc:
                    print("Discriminant was {:.3e}; set it to 0.".format(disc))
                newt0 = -b/2/a
            factor = self.ve2qf[(vtx.nr, el.nr)
                                ] if self.qf and dim > 1 else 1.0
            dt = (newt0 - t0) * factor
            dts.append(dt)
        if ret_dts:
            return dts
        else:
            return max(0, min(dts))

    # ------------------ #
    # diagnostic methods #
    # ------------------ #

    def checktent(self, tent, tau, tol=1.0e-15):
        """
        Verify that the norms of the gradients of all facets of the
        tent are less than the current wavespeed (and that the
        normal vectors used in computations are correct)
        """
        mesh = self.mesh
        gradnorms = []
        cbars = []
        vtx = tent.vertex
        for el in mesh[vtx].elements:
            gradtau = np.zeros(self.dim)
            cb = []
            for v in mesh[el].vertices:
                cb.append(self.cbar[v])
                gradtau += self.tau[v] * self.v2elgradphi[v.nr][el.nr]
            cbars.append(max(cb))
            gradnorms.append(la.norm(gradtau))

        if min(cbars) == max(cbars):
            # compare for the tent
            gradnorm = max(gradnorms)
            cbar = cbars[0]
            msg = "tent {} at {}, |∇ f|: {:.5f}, 1/c̅: {:.5f}, diff: {:.3e}{}"
            msg = msg.format(
                len(self.tents)-1, tent.vertex, gradnorm, 1./cbar,
                gradnorm - 1./cbar,
                ", at top" if tent.ttop == self.slabheight else "")
            if self.printcheck:
                print(msg)
            err = abs(gradnorm - 1./cbar)
            if err > self.maxerr and tent.ttop < self.slabheight:
                self.maxerr = err
            if gradnorm > self.maxgradtau:
                self.maxgradtau = gradnorm
            ignore = tent.ttop == self.slabheight or self.qf
            if not self.edge_based:
                assert err < tol or (ignore and gradnorm <= 1./cbar), \
                    "Max |∇ f| = 1/c̅ (or is less if at slab top). " + msg
        else:
            # need to compare separately for each facet
            for i, norm in enumerate(gradnorms):
                cbar = cbars[i]
                assert norm < 1./cbar + tol, \
                    "Max |∇ f| for facet <= 1/c̅ for facet"

            diffs = [abs(gn-1./cb) for (gn, cb) in zip(gradnorms, cbars)]
            md = min(diffs)
            ix = diffs.index(md)
            assert md < tol or tent.ttop == self.slabheight, \
                "There is a facet for which Max |∇ f| for tent equals 1/c̅" \
                + " (or is smaller if at slab top)"
            if self.printcheck:
                msg = "tent {} at {}, facet {}," \
                    + " |∇ f|: {:.5f}, 1/c̅: {:.5f}, diff: {:.3e}, {}"
                print(msg.format(
                    len(self.tents)-1, tent.vertex, ix, gradnorms[ix],
                    1./cbars[ix], gradnorms[ix] - 1./cbars[ix],
                    "at top" if tent.ttop == self.slabheight else ""))

    def ntentsforlevel(self):
        """
        Return the number of tents at each level
        """
        maxlevel = max(tent.level for tent in self.tents)
        results = []
        for lvl in range(maxlevel+1):
            results.append(
                (lvl, len([t for t in self.tents if t.level == lvl])))
        return results

    # ------------------------------ #
    # TentPitchedSlab Helper methods #
    # ------------------------------ #

    def tentpts(self, tent):
        """
        Get data for 1D Matplotlib
        """
        xs = [self.vpt[v][0] for v in self.v2v[tent.vertex]]
        ts = [tent.nbtime[i] for i in range(len(self.v2v[tent.vertex]))]
        vx = self.vpt[tent.vertex][0]

        if len(xs) == 1:
            xs = [vx, vx, xs[0]]
            ts = [tent.tbot, tent.ttop, ts[0]]
        else:
            xs = [xs[0], vx, xs[1], vx]
            ts = [ts[0], tent.tbot, ts[1], tent.ttop]
        return xs, ts

    def draw(self, n=None, step=False):
        """
        1D Matplotlib tent scene
        """
        def avg(lst):
            return sum(lst)/len(lst)
        colors = list('bgrcmyk')
        if step:
            plt.ion()
        if self.mesh.dim > 1:
            raise ValueError("Draw only supported for 1D")
        for i, t in enumerate(self.tents[:n]):
            xs, ts = self.tentpts(t)
            xpos, tpos = avg(xs), avg(ts)
            xs.append(xs[0])
            ts.append(ts[0])
            plt.plot(xs, ts, color=colors[t.level % 7])
            plt.text(xpos, tpos, "{}".format(i),
                     horizontalalignment='center',
                     verticalalignment='center')
            plt.ylim([0, self.slabheight*1.1])
            plt.xlim([0, 1])
            if step:
                input("enter")
        if not step:
            plt.show()

    def DrawPitchedTentsGL(self):
        """
        For 2D visualization using tentswebgui
        """
        nlevels = 0
        tentdata, tenttimes = [], []

        for i, tent in enumerate(self.tents):
            for el in tent.els:
                tentdata.extend([i, tent.level, tent.vertex.nr, el])
                if tent.level > nlevels:
                    nlevels = tent.level
                if self.mesh.dim == 2:
                    eid = ng.ElementId(ng.VOL, el)
                    vs = self.mesh[eid].vertices
                    for v in vs:
                        if v.nr in tent.nbv:
                            idx = tent.nbv.index(v.nr)
                            tenttimes.append(tent.nbtime[idx])
                        else:
                            tenttimes.append(tent.tbot)
                    tenttimes.append(tent.ttop)
        return tentdata, tenttimes, len(self.tents), nlevels+1

    def Draw3DTentPlt(self, tentnr):
        """
        Draw a single 3D tent using Matplotlib. Colors are used to represent
        the times of the neighbor vertices.  The tent pole height is represented
        by the size of the central vertex.
        """
        if self.dim != 3:
            raise NotImplementedError("Only supported for 3D spatial mesh")
        mesh = self.mesh
        tent = self.GetTent(tentnr)
        vtx = tent.vertex
        nbv = list(tent.nbv)
        tt = tent.ttop
        tb = tent.tbot
        nbtime = list(tent.nbtime)

        mvs = [vtx]+[ng.NodeId(ng.VERTEX, nb) for nb in nbv]
        mvs = [mesh[mv] for mv in mvs]
        print("mvs", mvs)
        mels = [mesh[ng.ElementId(ng.VOL, e)] for e in tent.els]
        pts = np.array([v.point for v in mvs])
        facetvs = [[mvs.index(v) for v in el.vertices if v != mvs[0]]
                   for el in mels]
        fig = plt.figure()
        fig.suptitle('Tent {} at level {}'.format(tentnr, tent.level))
        ax = fig.add_subplot(111, projection='3d')
        # edges from central vertex to neighbors
        for i in range(1, pts.shape[0]):
            ax.plot([pts[0, 0], pts[i, 0]], [pts[0, 1], pts[i, 1]],
                    [pts[0, 2], pts[i, 2]], color='blue', linewidth=.5)

        # outlines of facets
        for i, f in enumerate(facetvs):
            xs = [pts[v, 0] for v in f] + [pts[f[0], 0]]
            ys = [pts[v, 1] for v in f] + [pts[f[0], 1]]
            zs = [pts[v, 2] for v in f] + [pts[f[0], 2]]
            ax.plot(xs, ys, zs, color='red', linewidth=.5)

        # vertices
        s = [20 + 50*(tt-tb)] + [10]*(len(mvs)-1)
        scatter = ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2],
                             c=[tt]+nbtime, s=s)

        legend = ax.legend(*scatter.legend_elements(),
                           loc="lower left", title="Times")
        ax.add_artist(legend)
        print("Tent pole top: {:.3f}, bottom: {:.3f}, height: {:.3f}"
              .format(tb, tt, tt-tb))
        return ax


if __name__ == '__main__':

    # -------------------#
    # Simulation options #
    # -------------------#

    # slab height (in time)
    slabht = 0.5
    # spatial dimension
    dim = 2

    # initial ready_factor for greedy algorithms
    ready_factor = .5
    # if set, reset the ready_factor to its original value as soon as the
    # ready vertices collection has been repopulated
    reset_rf = False
    # if set, perform Greedy II - choose a vertex from ready_vs
    # with minimal layer and  maximal ktilde value
    maxready = False
    # if set, compute a mesh quality factor QF in the 2D case
    # to reduce computed kbar values
    qf = False
    # if set, use the edge_based method to compute dt
    edge_based = False
    # Using the edge-based method,
    # set Ctau to 1 for trial run and note maxgradtau
    # then set Ctau to speed/maxgradtau
    # for subsequent runs with the same mesh.
    Ctau = 1.0
    # draw tents in matplotlib in 1D, or generate HTML in 2D,
    # or print tents and draw single tent in matplotlib in 3D
    draw = True
    # tent number to draw for 3D case
    tentnr = 0

    # if set, use uniform (unstructured in 2D) mesh for 1D or 2D
    uniform = True
    # 1D number of elements for uniform mesh
    nelem1D = 20
    # if set, use torus in 3D case
    torus = False
    # number of refinements for 2D or 3D
    nrefinements = 2
    # save the generated mesh to a file
    save_mesh = False
    # name of saved mesh to use
    savedmesh = None
    # savedmesh = "d3_ref1.vol"
    # if set, draw mesh in netgen
    drawmesh = False

    # wavespeed and amplitude for coefficient functions
    speed = 1.0
    amplitude = 1.0
    # choice of wavespeed coefficient function (0=constant, 1=sine)
    cf = 0

    # Experimental increment to wavespeed applied each time a tent is created.
    #   The slab should complete without gradient errors in cases:
    #   accel<=0.00001 with slabheight 0.2, cf=0, speed=1, ready_factor=1/8
    #   and a twice refined unit_square(maxh=.2) mesh.
    accel = 0.0

    # if set, the max |∇ f| is computed over all facets of each new tent
    # and this value is compared with 1/c̅.  The difference should always
    # be zero except in the case when the tent top time equals slabheight,
    # in which case max |∇ f| < 1/c̅.
    chkgrad = True
    # tolerance for gradient check
    tol = 1.0e-10
    # if set, a line is printed for each tent showing the results of the
    # gradient check
    printcheck = False
    # if set, the discriminant is printed if it becomes negative
    printdisc = True

    # ---------------------- #
    # end simulation options #
    # ---------------------- #

    wavespd = None
    if dim == 1:
        if cf == 0:
            wavespd = ng.CoefficientFunction(speed)
        elif cf == 1:
            wavespd = speed*(amplitude*ng.sin(np.pi*ng.x) + 1)
    elif dim == 2:
        if cf == 0:
            wavespd = ng.CoefficientFunction(speed)
        elif cf == 1:
            wavespd = speed*(amplitude*(1-ng.sin(np.pi*ng.x)
                                        * ng.sin(np.pi*ng.y)) + 1)
    elif dim == 3:
        if cf == 0:
            wavespd = ng.CoefficientFunction(speed)
        elif cf == 1:
            wavespd = speed*(amplitude*ng.sin(np.pi*ng.x)
                             * ng.sin(np.pi*ng.y) * ng.sin(np.pi*ng.z) + 1)

    msh = None
    if savedmesh is None:
        if dim == 1:
            if uniform:
                msh = Make1DMesh([[0, 1]], [nelem1D], bc=[1, 3])
            else:
                pts = [0, 0.04, 0.08, 0.12, 0.18, 0.24, 0.30, 0.35,
                       0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
                       0.75, 0.80, 0.85, 0.90, 0.95, 1.0]
                msh = Make1DMeshSpecified(pts, bc=[1, 3])
        elif dim == 2:
            if uniform:
                from netgen.geom2d import unit_square
                msh = unit_square.GenerateMesh(maxh=0.2)
            else:
                from netgen.geom2d import SplineGeometry
                geo = SplineGeometry()
                geo.AddRectangle((0, 0), (1, 1))
                geo.AddCircle((0.5, 0.5), r=0.1, leftdomain=0, rightdomain=1)
                msh = geo.GenerateMesh(maxh=0.2)
        elif dim == 3:
            if torus:
                from netgen.meshing import Vec3d, Point3d
                from netgen.csg import CSGeometry, Torus
                geo = CSGeometry()
                tor = Torus(Point3d(0, 0, 0), Vec3d(0, 0, 1), 4, 1)
                geo.Add(tor)
                msh = geo.GenerateMesh(maxh=0.5)
            else:
                from netgen.csg import unit_cube
                msh = unit_cube.GenerateMesh(maxh=.3)

        msh = ng.Mesh(msh)
        if dim > 1:
            for i in range(nrefinements):
                msh.Refine()

        if dim == 3 and torus:
            msh.Curve(2)

        if save_mesh:
            msh.ngmesh.Save("d{}_ref{}.vol".format(dim, nrefinements))
    else:
        msh = ng.Mesh(savedmesh)

    if drawmesh:
        import netgen.gui
        ng.Draw(msh)
    slab = TentPitchedSlab(msh)
    tau, ktilde = slab.pitch_tents(slabht, wavespd, accel, chkgrad,
                                   tol, printcheck, printdisc,
                                   ready_factor, reset_rf, maxready, qf,
                                   edge_based, Ctau)
    if draw:
        if dim == 1:
            slab.draw()
        elif dim == 2:
            Draw(slab, "tents.html")
        else:
            # we don't have a good way to visualize a 3D tentslab yet.
            # but you can visualize a single tent
            slab.print_tents()
            ax = slab.Draw3DTentPlt(tentnr)
            plt.show()

    if min(tau.values()) < slabht:
        assert max(ktilde.values()) == 0
        msg = """
Tent slab is incomplete, but ktilde is zero for all vertices.
Perhaps the wavespeed is increasing too quickly or the spatial
variations in wavespeeed are too large for the specified slab height."""
        print(msg)
    print("max norm grad tau:", slab.maxgradtau)
    print("max error over all tents:", slab.maxerr)
    print("{} tents, {} layers".format(len(slab.tents),
                                       max(t.level for t in slab.tents)-1))
    print("Number of tents in each layer")
    print(slab.ntentsforlevel())
