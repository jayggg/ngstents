"""
Module erickson

A pure Python implementation the algorithm in the 2005 paper by
Erickson et. al. of a tent pitched slab
Some options are provided for experimentation and debugging

Visualizing the tent-pitched slab:
For 2D trials, running from the command line generates an HTML file.
For 3D trials, running in iPython allows visualization via pyplot
"""

from collections import defaultdict

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import ngsolve as ng

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
        self.cbar = None
        self.tau = None
        self.ready_vs = None
        self.tol = 1.0e-15
        self.printcheck = False
        self.printdisc = False
        self.maxerr = 0.0
        self.maxgradtau = 0.0
        self.v2oppfacetv = defaultdict(list)
        self.ve2trig = {}  # for 2D case
        self.ve2trigv = {}  # for 3D case
        self.eps1 = None
        self.eps2 = None

    def print_tents(self):
        for i, t in enumerate(self.tents):
            print("tent ", i, t, "\n")

    def GetTent(self, tentnr):
        return self.tents[tentnr]

    def set_opfacet_map(self):
        """
        construct a dict that maps each vertex to a list of its 'facets',
        where each 'facet' is represented by a list of vertices.
        """
        mesh = self.mesh
        for element in mesh.Elements():
            facets = set(mesh[f] for f in element.facets)
            for v in element.vertices:
                for f in facets:
                    if v not in f.vertices:
                        self.v2oppfacetv[v].append(list(f.vertices))

    def set_ve2trig(self):
        """
        construct a dict that maps a pair (central vertex, element) to
        a list of three vertices in which the central vertex appears first
        """
        for v in self.mesh.vertices:
            for e in self.mesh[v].elements:  # e is a trig
                mel = self.mesh[e]
                vs = [v] + [ev for ev in mel.vertices if ev != v]
                self.ve2trig[(v, e)] = vs

    def set_ve2trigv(self):
        """
        construct a dict that maps a pair (central vertex, element) to
        a list of pairs of the form (trig, other vtx)
        """
        for v in self.mesh.vertices:
            for e in self.mesh[v].elements:  # e is a tet
                mel = self.mesh[e]
                vs = [ev for ev in mel.vertices if ev != v]
                lst = [([v, vs[1], vs[2]], vs[0]),
                       ([v, vs[2], vs[0]], vs[1]),
                       ([v, vs[0], vs[1]], vs[2])]
                self.ve2trigv[(v, e)] = lst

    def pitch_tents(self, slabheight, wavespeed, eps1,
                    eps2, checkgrad=False, tol=1.0e-15,
                    printcheck=False, printdisc=False):
        """
        Generates all tents to fill the tent-pitched slab

        INPUTS

        slabheight: float - the time dimension of the slab to be generated
        wavespeed:  NGSolve CoefficientFunction - the initial wavespeed at
                      each spatial coordinate
        eps1: float in (0, 1/2] - the progress constraint constant used
                   by the nongreedy algorithm.  Set to None for greedy.
        eps2: float in (0, 1] - 3D case progress constraint constant used
                   by the nongreedy algorithm.  Set to None for greedy.
        checkgrad: optional boolean - print the computed grad norm, 1/cbar
                   and the error for each facet based on the dt value just
                   computed.
        tol: float (default 1.0e-14) - tolerance for gradient check.
        printcheck: optional boolean - if set, prints a line in gradient check
        """
        self.slabheight = slabheight
        self.wavespeed = wavespeed
        self.eps1 = eps1
        self.eps2 = eps2
        self.checkgrad = checkgrad
        self.tol = tol
        self.printcheck = printcheck
        self.printdisc = printdisc

        self.set_opfacet_map()
        if self.dim == 2:
            self.set_ve2trig()
        else:
            self.set_ve2trigv()

        # initial advancing front at vertices
        tau = self.tau = {v: 0 for v in self.vertices}
        # wavespeed for each vertex
        self.cbar = {v: wavespeed(self.mesh(*self.vpt[v]))
                     for v in self.vertices}
        if min(self.cbar.values()) <= 1e-16:
            raise ValueError("Wavespeed is not strictly positive")

        # initial level for each vertex
        vs_level = {v: 0 for v in self.vertices}

        # initially every vertex is a local minimum
        self.ready_vs = self.vertices.copy()
        mintau = 0

        while mintau < slabheight:
            # select a tent vertex and remove it from ready_vs
            minlevel = min(vs_level[v] for v in self.ready_vs)
            levelvs = [v for v in self.ready_vs if vs_level[v] == minlevel]
            vi = levelvs[0]
            self.ready_vs.remove(vi)

            # pitch tent, update tau, ready_vs, ktilde if greedy method
            t = self.make_tent(vi, slabheight, vs_level)
            if t is None:
                return tau

            if checkgrad:
                self.checktent(t, tau, tol=self.tol)
                if self.dim == 3:
                    self.checkgradHtau(t, tau, tol=self.tol)

            mintau = min(tau.values())
        return tau

    def make_tent(self, vi, slabheight, vs_level):
        """
        Generate a tent using the Erickson algorithm
        and update ready_vs.
        """
        tau = self.tau
        dt = self.compute_dt(vi)
        if dt <= 0:
            print("len(self.tents)", len(self.tents))
            print("self.vpt[vi]", self.vpt[vi])
            print("minimal vertex should have positive dt")
            # raise ValueError("minimal vertex should have positive dt")
            return None

        # construct the tent and append it
        t = Tent(vi, tau[vi],
                 min(slabheight, tau[vi]+dt), vs_level[vi],
                 [e.nr for e in self.mesh[vi].elements])
        self.tents.append(t)

        # update the advancing front
        tau[vi] = t.ttop

        # update the level for the central vertex
        vs_level[vi] += 1
        for nb in self.v2v[vi]:
            t.nbv.append(nb.nr)
            t.nbtime.append(tau[nb])
            # update the levels for neighbor vertices
            if vs_level[nb] < t.level + 1:
                vs_level[nb] = t.level + 1

        # add any neighbors which are now local mins to ready_vs
        for nb in self.v2v[vi]:
            if tau[nb] < self.slabheight and nb not in self.ready_vs and \
                    tau[nb] <= min(tau[vv] for vv in self.v2v[nb]):
                self.ready_vs.append(nb)
        return t

    def compute_dt(self, vtx):
        eps2 = self.eps2
        dts = []
        if self.dim == 2:
            for e in self.mesh[vtx].elements:
                trig = self.ve2trig[(vtx, e)]
                # make spacetime points
                pts = [list(self.vpt[v])+[self.tau[v]] for v in trig]
                cbar = max((self.cbar[v] for v in trig))
                dts.append(self.get_dt(pts, cbar))
        elif self.dim == 3:
            for i, e in enumerate(self.mesh[vtx].elements):
                # 2D cone constraints and progress constraint for each edge
                # containing our vertex
                lst = self.ve2trigv[(vtx, e)]
                for j, item in enumerate(lst):
                    trig, po = item
                    pts = [list(self.vpt[v])+[self.tau[v]] for v in trig]
                    po = list(self.vpt[po])
                    cbar = max((self.cbar[v] for v in trig))
                    dt = self.get_dt(pts, cbar, po=po, vtx=vtx, e=e, trigix=j)
                    dts.append(dt)
            # Now ensure that the 3D cone constraint is not exceeded
            for i, facet in enumerate(self.v2oppfacetv[vtx]):
                # make spacetime points
                pts = np.array([list(self.vpt[v])+[self.tau[v]]
                                for v in facet])
                cbar = max((self.cbar[v] for v in facet))
                t0 = self.tau[vtx]
                p = self.vpt[vtx]
                pH = calc_pH(pts[:, :-1], p)
                pHt = calc_pHt(pts, pH)
                sigmaF = calc_sigmaF(pts[:, :-1], p)
                term = eps2 + (1-eps2)*np.sqrt(1-sigmaF**2)
                newt0 = pHt + la.norm(p-pH)*term/cbar
                dts.append(newt0 - t0)
        else:
            raise ValueError(
                "Only spatial dimensions 2 and 3 are supported")
        return max(0, min(dts))

    def get_dt(self, pts, cbar, po=None, vtx=None, e=None, trigix=None):
        """
        INPUTS

        pts: the spacetime coordinates of the central vertex followed by
             those of the vertices of an opposite facet (edge)
             the spatial coordinates may be in 2D or 3D
        cbar: the wavespeed for the current element
        po: an optional vertex not in the plane of the other vertices

        Returns the optimal time increment for the central vertex
        """
        A = np.array(pts)

        tp, tq, tr = A[:, -1]
        if tr < tq:
            A = A[[0, 2, 1], :]
            tp, tq, tr = A[:, -1]

        p, q, r = A[:, :-1]
        pq, rq, trq = p-q, r-q, tr-tq
        nrq = la.norm(rq)
        wp = la.norm(np.cross(pq, rq)) / nrq
        if po is None:
            disc = (nrq/cbar)**2 - trq**2
        else:
            sigmaF = calc_sigmaF(A[:, :-1], po)
            # set dt1 so that ∥∇_H tau∥ = (1-eps2)sigmaF/cbar (checked)
            disc = (nrq*sigmaF*(1-self.eps2)/cbar)**2 - trq**2

        if disc < 0:
            print("negative disc", disc, "setting to zero")
            disc = 0
        conetp = tq + trq*pq.dot(rq)/nrq**2 + wp*np.sqrt(disc)/nrq

        if (False and self.dim == 3 and len(self.tents) == 86
                and vtx.nr == 89 and e.nr == 672 and trigix == 0):
            gradHtaunorm = calc_gradHtaunorm(A)
            pH = calc_pH(A[:, :-1], po)
            pF = closest_pt_to_Trig(A[:, :-1], po)
            constraint = (1-self.eps2)*sigmaF/cbar
            npq = la.norm(pq)
            pr = p - r
            npr = la.norm(pr)
            theta = np.arccos(pq.dot(rq)/npq/nrq)
            alpha = np.arcsin(wp/npr)
            gradHtaunorm2 = tr/nrq/np.sin(theta)
            check = pq.dot(rq)/npq*sigmaF*(1-self.eps2)
            print("check", check)
            print("A", A)
            print("vtx", vtx)
            print("e", e)
            print("trigix", trigix)
            print("gradhtaunorm before", gradHtaunorm)
            print("gradhtaunorm2 before", gradHtaunorm2)
            print("constraint", constraint)
            print("tp", tp)
            print("tq", tq)
            print("tr", tr)
            print("p", p)
            print("q", q)
            print("r", r)
            print("po", po)
            print("sigmaF", sigmaF)
            print("pH", pH)
            print("pF", pF)
            print("trq", trq)
            print("disc", disc)
            print("pq.dot(rq)", pq.dot(rq))
            print("nrq", nrq)
            print("npq", npq)
            print("nrq**2", nrq**2)
            print("wp", wp)
            print("trq*pq.dot(rq)/nrq**2", trq*pq.dot(rq)/nrq**2)
            print("cone1", conetp)
            print("theta", theta)
            print("alpha", alpha)

        prog_constr = tr + (1-self.eps1)*wp/cbar
        newtp = min(conetp, prog_constr)

        dtp = newtp - tp
        return dtp

    # ------------------ #
    # diagnostic methods #
    # ------------------ #

    def export_tau(self, ktilde, vert, dts):
        """
        temporarily provide export to npz for testing
        """
        ar = np.array([[v.nr] + list(self.vpt[v]) + [time] + [ktilde[v]]
                       for v, time in self.tau.items()])
        v2v = {k.nr: [v.nr for v in lst] for k, lst in self.v2v.items()}
        facets = [[v.nr for v in facet] for facet in self.v2oppfacetv[vert]]
        np.savez('taupts', ar=ar, v2v=v2v, vert=vert.nr, dts=dts,
                 facets=facets)

    def checkgradHtau(self, tent, tau, tol=1.0e-15):
        """
        Verify that the norm of the gradient of the restriction of tau
        to F is less than or equal to (1-eps2)sigmaF for each face F of
        the tent that includes the tent vertex.
        """
        ok = True
        vtx = tent.vertex
        for i, e in enumerate(self.mesh[vtx].elements):
            # 2D cone constraints and progress constraint for each edge
            # containing our vertex
            lst = self.ve2trigv[(vtx, e)]
            for j, item in enumerate(lst):
                trig, po = item
                pts = np.array([list(self.vpt[v])+[self.tau[v]] for v in trig])
                po = np.array(list(self.vpt[po]))
                # make spacetime points
                cbar = max((self.cbar[v] for v in trig))
                sigmaF = calc_sigmaF(pts[:, :-1], po)
                gradHtau = calc_gradHtaunorm(pts)
                if gradHtau > (1-self.eps2)*sigmaF/cbar + tol:
                    print("progress constraint violated")
                    print("gradHtau", gradHtau)
                    print("(1-self.eps2)*sigmaF/cbar + tol",
                          (1-self.eps2)*sigmaF + tol)
                    print("po", po)
                    print("sigmaF", sigmaF)
                    print("self.eps2", self.eps2)
                    print("len(self.tents)", len(self.tents))
                    print("element", i)
                    print("trig", j)
                    print("pts", pts)
                    ok = False
        return ok

    def checktent(self, tent, tau, tol=1.0e-15):
        """
        Verify that the norms of the gradients of all facets of the
        tent are less than the current wavespeed (and that the
        normal vectors used in computations are correct)
        """
        gradnorms = []
        cbars = []
        # a facet is represented by a list of its vertices
        facets = self.v2oppfacetv[tent.vertex]
        for i, facet in enumerate(facets):
            # make spacetime points
            verts = [tent.vertex, *facet]
            pts = [list(self.vpt[v])+[self.tau[v]] for v in verts]
            # get cbars and gradnorms
            cbars.append(max((self.cbar[v] for v in verts)))
            if self.mesh.dim == 2:
                gradnorms.append(compute_normgrad2D(pts))
            elif self.mesh.dim == 3:
                gradnorms.append(compute_normgrad3D(pts))

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
            if err > self.maxerr and tent.ttop != self.slabheight:
                self.maxerr = err
            if gradnorm > self.maxgradtau:
                self.maxgradtau = gradnorm
            if gradnorm >= 1./cbar + tol:
                print("gradnorm", gradnorm)
                print("# of tents", len(self.tents))
                print("gradnorms", gradnorms)

            assert gradnorm < 1./cbar + tol
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

    def ntentsforlevel(self):
        maxlevel = max(tent.level for tent in self.tents)
        results = []
        for lvl in range(maxlevel+1):
            results.append(
                (lvl, len([t for t in self.tents if t.level == lvl])))
        return results

# --------------#
# Helper Methods
# --------------#


def calc_gradHtau(A):
    """
    INPUTS

    A: 3x4 array of spacetime 3D points

    returns ∇H τ, the gradient of the restriction of τ
    to the plane of the triangle in the form of coefficients μ and ν
    and unit vectors v̅ and n̅ along with t(q)
    """
    V = (A - A[0, :])[1:, :]
    tp, tr = V[:, -1]  # tp = t(p)-t(q), tr = t(r)-t(q)
    p, r = V[:, :-1]   # p̅ = p - q, r̅ = r - q
    nmr = la.norm(V[1, :-1])
    nmr2 = nmr*nmr
    pdotr = p.dot(r)
    wp = la.norm(np.cross(p, r))/nmr  # distance from p to opposite edge
    v = r/nmr
    w = pdotr*r/nmr2
    n = (p - w)/la.norm(p - w)
    mu = tr/nmr
    nu = (nmr2*tp - pdotr*tr)/wp/nmr2
    return mu, nu, v, n


def calc_pHt(A, pH):
    """
    INPUTS

    A: 3x4 array of spacetime 3D points representing the opposite face to p
    pH: The projection of p onto the plane of its opposite face

    returns tau(pH)
    """
    tq = A[0, -1]
    s = pH - A[0, :-1]
    mu, nu, v, n = calc_gradHtau(A)
    return tq + mu*v.dot(s) + nu*n.dot(s)


def calc_gradHtaunorm(A):
    """
    INPUTS

    A: 3x4 array of spacetime 3D points

    returns ∥∇H τ∥ the norm of the gradient of the restriction of τ
    to the plane of the triangle
    """
    mu, nu, _, _ = calc_gradHtau(A)
    return np.sqrt(mu*mu + nu*nu)


def calc_sigmaF(A, p):
    """
    INPUTS

    A: 3x3 array of spatial 3D points representing a triangle
    p: a 3D point not in the plane of the triangle
    """
    pH = calc_pH(A, p)
    pF = closest_pt_to_Trig(A, pH)
    sigmaF = la.norm(p-pH)/la.norm(p-pF)

    assert sigmaF <= 1 + 1.0e-14
    return min(sigmaF, 1.)


def calc_pH(A, p):
    """
    Compute the projection of p onto the plane of its opposite face

    INPUTS

    A: 3x3 array of spatial 3D points representing a triangle
    p: a 3D point not in the plane of the triangle
    """
    p = np.array(p)
    v0 = A[0, :]
    w1, w2 = (A - v0)[1:, :]
    n = np.cross(w1, w2)
    n = n/la.norm(n)
    pH = p - n.dot(p-v0)*n

    return pH


def closest_pt_to_Trig(A, p):
    """
    INPUTS

    pH: Projection of p to the plane containing the triangle F
    v0, v1, v2: vertices of triangle F
    w0, w1, w2: vectors def by w0 = v1-v0, w1 = v2-v1, w2 = v0-v2
    w0sq, w1sq, w2sq: squared norms of vectors
    """
    p = np.array(p)
    v0, v1, v2 = A

    w10, w21, w02 = np.roll(A, -1, 0) - A
    w20, w01, w12 = -w02, -w10, -w21

    w01sq, w12sq, w20sq = w01.dot(w01), w12.dot(w12), w20.dot(w20)
    # compute projection pH of p onto plane of triangle with verts in A
    pH = calc_pH(A, p)

    p0, c0 = closest_pt(pH, v0, v2, w10, w20, w20sq)
    p1, c1 = closest_pt(pH, v1, v0, w21, w01, w01sq)
    p2, c2 = closest_pt(pH, v2, v1, w02, w12, w12sq)
    if p0 is None and p1 is None and p2 is None:
        # inside triangle
        return pH
    if p0 is not None and 0 <= c0 <= 1:
        # on edge w20
        return p0
    if p1 is not None and 0 <= c1 <= 1:
        # on edge w01
        return p1
    if p2 is not None and 0 <= c2 <= 1:
        # on edge w12
        return p2
    if p0 is not None:
        # v0 if c0 < 0 else v2
        return p0
    if p1 is not None:
        # v1 if c1 < 0 else v0
        return p1
    # v2 if c2 < 0 else v1
    return p2


def closest_pt(p, v0, v2, w10, w20, w20sq):
    """
    INPUTS

    p, v1, v2: 3D points in same plane
    w10 = v1-v0
    w20 = v2-v0
    w20sq = norm squared of w20
    """
    # check if p and v1 on same side of w20
    d = np.cross(p-v0, w20).dot(np.cross(w10, w20))
    if d <= 0:
        # opposite side so check projection
        c = (p-v0).dot(w20)/w20sq
        # print("c", c, "w20sq", w20sq)
        if c < 0:
            # print("returning v0")
            return v0, c
        elif c > 1:
            # print("returning v2")
            return v2, c
        # print("returning projection")
        return v0 + c*w20, c  # projection
    else:
        return None, None


def compute_normgrad2D(pts, nonorm=False):
    """
    Given a list of three spacetime points, compute the normal
    to the plane they lie in and its gradient
    """
    A = np.array(pts)
    m20, m21, m22, V = A2m2D(A)
    a = la.det(m20)
    b = - la.det(m21)
    c = la.det(m22)
    if nonorm:
        return np.array([-a/c, -b/c])
    normgrad = np.sqrt((a*a + b*b)/(c*c))
    checknormal = V@np.array([a, b, c])
    nmlnorm = la.norm(checknormal)
    assert nmlnorm < 1.0e-15, "normal is orthogonal to plane"
    return normgrad


def A2m2D(A):
    """
    Get matrix of point differences and minor matrices
    """
    V = A[1:, :] - A[:-1, :]
    return V[:, [1, 2]], V[:, [0, 2]], V[:, [0, 1]], V


def compute_normgrad3D(pts):
    """
    Given a list of four spacetime points, compute the normal
    to the hyperplane plane they lie in and its gradient
    """
    A = np.array(pts)
    m30, m31, m32, m33, V = A2m3D(A)
    a, b, c, e = -la.det(m30), la.det(m31), -la.det(m32), la.det(m33)
    normgrad = np.sqrt((a*a + b*b + c*c)/(e*e))
    checknormal = V@np.array([a, b, c, e])
    nmlnorm = la.norm(checknormal)
    assert nmlnorm < 1.0e-15, "normal is orthogonal to hyperplane"
    return normgrad


def A2m3D(A):
    """
    Get matrix of point differences and minor matices
    """
    V = A[1:, :] - A[:-1, :]
    return V[:, [1, 2, 3]], V[:, [0, 2, 3]], \
        V[:, [0, 1, 3]], V[:, [0, 1, 2]], V


if __name__ == '__main__':

    # -------------------#
    # Simulation options #
    # -------------------#

    # slab height (in time)
    slabht = 0.5
    # spatial dimension 2 or 3
    dim = 2

    # parameter for progress constraint
    eps1 = 0.1
    # parameter for 3D constraint
    eps2 = 0.3  # if eps2 >= 0.46 the algorithm gets stuck on tent 86
    # generate HTML in 2D,
    # or print tents and draw single tent in matplotlib in 3D
    draw = False
    # tent number to draw for 3D case
    tentnr = 0

    # if set, use uniform unstructured mesh for 2D
    uniform = True
    # if set, use torus in 3D case
    torus = False
    # number of refinements for 2D or 3D
    nrefinements = 2
    # save the generated mesh to a file
    save_mesh = False
    # name of saved mesh to use
    savedmesh = None
    # savedmesh = "d2_ref2.vol"
    # if set, draw mesh in netgen
    drawmesh = False

    # wavespeed and amplitude for coefficient functions
    speed = 1.0
    amplitude = 1.0
    # choice of wavespeed coefficient function (0=constant, 1=sine)
    cf = 0

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
        if dim == 2:
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
    tau = slab.pitch_tents(slabht, wavespd, eps1, eps2, chkgrad,
                           tol, printcheck, printdisc)
    if draw:
        if dim == 2:
            Draw(slab, "tents.html")
        else:
            # we don't have a good way to visualize a 3D tentslab yet.
            # but you can visualize a single tent
            slab.print_tents()
            ax = slab.Draw3DTentPlt(tentnr)
            plt.show()

    if min(tau.values()) < slabht:
        msg = """ Tent slab is incomplete -- why?"""
        print(msg)
    print("max norm grad tau:", slab.maxgradtau)
    print("max error over all tents:", slab.maxerr)
    print("{} tents, {} layers".format(len(slab.tents),
                                       max(t.level for t in slab.tents)-1))
    print("Number of tents in each layer")
    print(slab.ntentsforlevel())
