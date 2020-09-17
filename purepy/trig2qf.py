"""
A module for testing Quality Factor variants of the Greedy I algorithm
"""
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


class simplex(object):
    """
    A simplex
    """

    def __init__(self, pts, cbar=1.0, qalg=2):
        """
        INPUTS

        pts: d+1 x d numpy array, each row a vertex of a conforming simplex
        cbar: float wavespeed
        """
        self.pts = pts
        self.d = pts.shape[1]
        self.cbar = cbar
        self.qalg = qalg
        self.nmls = None
        self.ts = np.zeros(self.d+1)
        Ms, Msums, Mdot1s, Ainvs = [], [], [], []
        for i in range(self.d+1):
            At = np.roll(pts - pts[i, :], -(i+1), 0)[:self.d, :]
            Ainvs.append(la.inv(At.T))
            M = la.inv(At@At.T)
            Ms.append(M)
            Mdot1s.append(M.sum(axis=0))
            Msums.append(M.sum())
        self.Ms, self.Msums, self.Mdot1s = Ms, Msums, Mdot1s
        self.Ainvs = Ainvs
        self.edgelens = None
        self.q = None

    def get_edgelens(self):
        """
        lengths of edges of the trig
        """
        if self.edgelens is None:
            v = np.roll(self.pts, -1, 0) - self.pts
            self.edgelens = np.roll((((v)**2).sum(1))**.5, -1)
        return self.edgelens

    def get_qv1(self):
        """
        mesh quality factor version 1
        """
        if self.q is None:
            edgelens = self.get_edgelens()
            self.q = [min(edgelens)/max(edgelens)]*3
        return self.q

    def get_qv2(self):
        """
        mesh quality factor version 2
        """
        if self.q is None:
            self.q = []
            edgelens = self.get_edgelens()
            for i in range(3):
                j, k = (i+1) % 3, (i+2) % 3
                self.q.append(
                    min(1., edgelens[i]/max(edgelens[j], edgelens[k])))
        return self.q

    def gtdotn(self, ix):
        """
        Check if the condition ∇ τ ⋅ nᵢₓ > 0
        holds for the vertex ix

        INPUTS

        ix: vertex index to check
        """
        return self.gradtau().dot(self.normals()[ix])

    def normals(self):
        """
        Compute the normal vectors to each facet indexed by opposite vertex
        Caching the computed result
        """
        if self.nmls is None:
            # normal vector of facet opposite v0 for unit tet
            nml0 = np.ones(self.d)
            nmls = []
            for i in range(self.d + 1):
                nmls.append(nml0.dot(self.Ainvs[i]))
            self.nmls = np.array(nmls)
        return self.nmls

    def gradtau(self):
        """
        Compute ∇ τ for the given vertex times

        INPUTS

        ts: numpy array with d+1 elements representing vertex times
        """
        dts = self.dts(0)
        return dts.dot(self.Ainvs[0])

    def gradnormsq(self):
        """
        Compute ||∇ τ||²

        INPUTS

        ts: numpy array with d+1 elements representing vertex times
        """
        dts = self.dts(0)
        return dts.dot(self.Ms[0]).dot(dts)

    def dts(self, ix):
        """
        Compute dtᵢₓ, the vector of time differences relative to vᵢₓ

        INPUTS

        ts: numpy array with d+1 elements representing vertex times
        ix: index of reference vertex
        """
        return np.roll(self.ts - self.ts[ix], -(ix+1))[:self.d]

    def kbar(self, ix):
        """
        Compute k̅ᵢₓ for vertex ix

        INPUTS

        ts: numpy array with d+1 elements representing vertex times
        """
        if self.qalg == 1:
            q = self.get_qv1()  # quality factor version 1
        elif self.qalg == 2:
            q = self.get_qv2()  # quality factor version 2
        else:
            q = [1.]*3          # no quality factor (Greedy I)
        M, Msum, Mdot1 = self.Ms[ix], self.Msums[ix], self.Mdot1s[ix]
        c = self.cbar
        dt = self.dts(ix)
        m1dt = Mdot1.dot(dt)
        disc = m1dt**2 - Msum * (dt.dot(M).dot(dt) - 1/c**2)

        if disc < 0:
            disc = 0.
        return (m1dt + np.sqrt(disc))/Msum*q[ix]


class slab2d(object):
    """
    A tent slab based on a 2D simplicial mesh with 2 triangles
    """

    def __init__(self, cbar=1.0, vpts=None, qalg=2, maxktilde=False,
                 edge_based=False, Ctau=1.0):
        """
        INPUTS
        cbar: float wavespeed default 1.0
        vpts: 4x2 numpy array representing 4 2D points (if not set
              (if not specified, random vertices are generated)
        qalg: 0: no QF algorithm, 1: first QF algorithm, 2: second QF algorithm
        maxktilde: If true, use Greedy II, else Greedy I
        """
        self.cbar = cbar
        self.qalg = qalg
        self.maxktilde = maxktilde
        self.edge_based = edge_based
        self.Ctau = Ctau
        self.edgetimes = None
        if vpts is None:
            pts1 = np.random.rand(3, 2)
            diffs1 = (pts1 - pts1[0, :])[1:, :]
            if la.det(diffs1) <= 0:
                pts1[[0, 1], :] = pts1[[1, 0], :]

            pts2i = pts1[[0, 2], :]
            pts2 = np.vstack((pts2i, np.random.rand(2)))
            diffs2 = (pts2 - pts2[0, :])[1:, :]
            while la.det(diffs2) <= 0:
                pts2 = np.vstack((pts2i, np.random.rand(2)))
                diffs2 = (pts2 - pts2[0, :])[1:, :]
            self.vpts = np.vstack((pts1, pts2[-1, :]))
        else:
            self.vpts = vpts
            pts1 = self.vpts[:-1, :]
            pts2 = self.vpts[[0, 2, 3], :]

        self.v2el = {0: [0, 1], 1: [0], 2: [0, 1], 3: [1]}
        self.el2v = {0: [0, 1, 2], 1: [0, 2, 3]}
        self.v2v = {0: [1, 3], 1: [0, 2], 2: [1, 3], 3: [0, 2]}
        self.edges = [0, 1, 2, 3, 4]
        self.v2edges = {0: [0, 2, 4], 1: [0, 1], 2: [1, 2, 3], 3: [3, 4]}
        self.edge2vs = {0: [0, 1], 1: [1, 2], 2: [0, 2], 3: [2, 3], 4: [0, 3]}
        self.els = [simplex(pts1, cbar, qalg), simplex(pts2, cbar, qalg)]
        self.ts = np.zeros(4)
        self.minstep = 1.0e16
        self.levels = [0, 0, 0, 0]
        self.vref = self.ktildes()
        self.ready_factor = 0.5
        self.maxgradtau = 0.0

    def elix(self, e, v):
        """
        Return the index of v in the context of element e
        """
        vs = self.el2v[e]
        if v not in vs:
            raise ValueError()
        return vs.index(v)

    def updatets(self, v):
        """
        Update the time for the specified vertex based on global times
        """
        for e in self.v2el[v]:
            ix = self.elix(e, v)
            self.els[e].ts[ix] = self.ts[v]

    def maxgradnormsq(self):
        return max(el.gradnormsq() for el in self.els)

    def gtdotns(self, v):
        return [np.round(self.els[e].gtdotn(self.elix(e, v)), 3)
                for e in self.v2el[v]]

    def allgtdotns(self):
        return [self.gtdotns(v) for v in range(4)]

    def isoptimal(self, vtx):
        """
        a vertex is optimal if the
        gtdotn corresponding to its minimal kbar is positive.
        """
        gtdotns = self.gtdotns(vtx)
        ix = np.argmin(self.kbars(vtx))
        return gtdotns[ix] > 0

    def kbars(self, v):
        """
        Given a vertex, compute time increments which would set ||∇ τ|| = 1/cbar
        for each adjacent element.
        """
        result = []
        for e in self.v2el[v]:
            el = self.els[e]
            eix = self.elix(e, v)
            result.append(el.kbar(eix))
        return result

    def ktildes(self):
        """
        Compute the time increment for each vertex v
        """
        result = []
        for v in range(4):
            if self.edge_based:
                # min of vertex time difference + edgetime
                vt = self.ts[v]
                dts = []
                for e in self.v2edges[v]:
                    v1 = [vv for vv in self.edge2vs[e] if vv != v][0]
                    dts.append(self.ts[v1]-vt + self.get_edgetimes()[e])
                result.append(min(dts))
            else:
                # min over elements incident on v of ||∇ τ|| = 1/cbar
                kbars = self.kbars(v)
                result.append(max(0, min(kbars)))
        return result

    def get_edgetimes(self):
        """
        Compute the edge times for the edge-based algorithm
        """
        if self.edgetimes is None:
            self.edgetimes = []
            for e in self.edges:
                v1, v2 = self.edge2vs[e]
                length = (((self.vpts[v1] - self.vpts[v2])**2).sum())**0.5
                self.edgetimes.append(length/self.cbar * self.Ctau)
        return self.edgetimes

    def step(self, verbose=False):
        """
        Implements Greedy I and Greedy II
        Select a vertex with minimal layer from readyvs
        Pitch a tent at that vertex
        """
        allktildes = self.ktildes()
        if max(allktildes) < 1.e-15:
            # print("all ktildes are zero")
            return max(allktildes)
        readyvs = [i for i in range(4)
                   if allktildes[i] >= self.ready_factor*self.vref[i]]
        while max(allktildes) > 1.e-16 and len(readyvs) == 0:
            print("reduced ready factor")
            self.ready_factor /= 2.
            readyvs = [i for i in range(4)
                       if allktildes[i] >= self.ready_factor*self.vref[i]]
        minlevel = min(self.levels[v] for v in readyvs)
        lvlktildes = [kt for i, kt in enumerate(allktildes)
                      if self.levels[i] == minlevel and i in readyvs]
        if self.maxktilde:
            # Greedy II
            ktilde = max(lvlktildes)
        else:
            # Greedy I
            ktilde = lvlktildes[0]
        ix = allktildes.index(ktilde)
        levelsbefore = self.levels.copy()
        self.levels[ix] += 1
        for v in self.v2v[ix]:
            self.levels[v] = self.levels[ix]

        if verbose:
            print("vertex times before", self.ts)
            print("levels before", levelsbefore)
            print("ktildes before", allktildes)
            print("readyvs", readyvs)

        if ktilde < self.minstep:
            self.minstep = ktilde
        self.ts[ix] += ktilde
        self.updatets(ix)
        if self.maxgradnormsq() > self.maxgradtau:
            self.maxgradtau = self.maxgradnormsq()
        if verbose:
            print("updated vertex", ix)
            print("vertex times after", self.ts)
            print("levels after", self.levels)
            print("minstep", self.minstep)
            # preview ktildes for next step
            # ktildes = self.ktildes()
            # print("next step ktildes", self.ktildes)
        return ktilde

    def plot(self):
        """
        Draw the 2 triangle mesh
        """
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        for el in self.els:
            lpts = np.vstack((el.pts, el.pts[0, :])).tolist()
            ppts = list(zip(*lpts))
            ax.plot(ppts[0], ppts[1], 'g-')
        vpts = self.vpts.tolist()
        ppts = list(zip(*vpts))
        ax.scatter(ppts[0], ppts[1], c='blue')
        for i, pt in enumerate(vpts):
            plt.annotate("  v{}".format(i), pt)
        plt.show()


def rand_simplex(d=3, cbar=1.0):
    """
    Generate a random d-simplex with vertices ordered so that it is
    homeomorphic to the unit d-simplex
    """
    pts = np.random.rand(d+1, d)
    vs = np.roll(pts - pts[0, :], -1, axis=0)[:d, :]
    if la.det(vs) < 0:
        print("reordered to make homeomorphic to unit d-simplex")
        pts[[0, 1], :] = pts[[1, 0], :]
    return simplex(pts, cbar)


def unit_simplex(d=3, cbar=1.0):
    """
    Construct the unit d-simplex with given wavespeed
    """
    pts = np.zeros((d+1, d))
    pts[1:, :] = np.eye(d)
    return simplex(pts, cbar)


def gen(verbose=False, pts=None, qalg=2, maxktilde=False,
        edge_based=False, Ctau=1.0):
    if pts is None:
        s = slab2d(qalg=qalg, maxktilde=maxktilde, edge_based=edge_based,
                   Ctau=Ctau)
    else:
        s = slab2d(vpts=pts, qalg=qalg, maxktilde=maxktilde,
                   edge_based=edge_based, Ctau=Ctau)
    for i in range(1000):
        print()
        print("tent", i)
        ktmax = s.step(verbose=verbose)
        if ktmax < 1.0e-15:
            print("ktilde is zero; tent is stuck")
            break
    print("times", s.ts)
    return s, i


def multigen(qalg=0, maxktilde=False, edge_based=False, Ctau=1.0):
    while True:
        print()
        print("new mesh")
        s, nsteps = gen(verbose=False, qalg=qalg, maxktilde=maxktilde,
                        edge_based=edge_based, Ctau=Ctau)
        if nsteps < 10:
            break
    s.plot()
    return s, nsteps


if __name__ == "__main__":
    s, nsteps = multigen(maxktilde=True, edge_based=True)

    # qalg 2 does not ever seem to get stuck since we implemented the full
    # features of the Greedy II algorithm.

    # Example for final report
    # s,i = gen(True, pts=pts, qalg=0, maxktilde=False)
    # pts = np.array([[0.41343681, 0.15895487],
    #                 [0.71252672, 0.97572707],
    #                 [0.53145611, 0.72308145],
    #                 [0.05134836, 0.05218915]])
    # pts = np.array([[0.413437, 0.158955],
    #                 [0.712527, 0.975727],
    #                 [0.531456, 0.723081],
    #                 [0.051348, 0.052189]])
    # pts = np.array([[0.4134, 0.1590],
    #                 [0.7125, 0.9757],
    #                 [0.5315, 0.7231],
    #                 [0.0513, 0.0522]])

    # s, i = gen(True, pts=pts, qalg=0, maxktilde=False, edge_based=True)
