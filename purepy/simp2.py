import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


class simplex(object):
    """
    A simplex
    """

    def __init__(self, pts, cbar=1.0):
        """
        INPUTS

        pts: d+1 x d numpy array, each row a vertex of a conforming simplex
        cbar: float wavespeed
        """
        self.pts = pts
        self.d = pts.shape[1]
        self.cbar = cbar
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
        M, Msum, Mdot1 = self.Ms[ix], self.Msums[ix], self.Mdot1s[ix]
        c = self.cbar
        dt = self.dts(ix)
        m1dt = Mdot1.dot(dt)
        disc = m1dt**2 - Msum * (dt.dot(M).dot(dt) - 1/c)
        if disc < 0:
            disc = 0.
        return (m1dt + np.sqrt(disc))/Msum


class slab2d(object):
    """
    A tent slab based on a 2D simplicial mesh with 2 triangles
    """

    def __init__(self, cbar=1.0, vpts=None):
        """
        INPUTS
        cbar: float wavespeed default 1.0
        vpts: 4x2 numpy array representing 4 2D points (if not set
              (if not specified, random vertices are generated)
        """
        self.cbar = cbar
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

        self.v2e = {0: [0, 1], 1: [0], 2: [0, 1], 3: [1]}
        self.e2v = {0: [0, 1, 2], 1: [0, 2, 3]}

        self.els = [simplex(pts1, cbar), simplex(pts2, cbar)]
        self.ts = np.zeros(4)
        self.minstep = 1.0e16

    def elix(self, e, v):
        """
        Return the index of v in the context of element e
        """
        vs = self.e2v[e]
        if v not in vs:
            raise ValueError()
        return vs.index(v)

    def updatets(self, v):
        """
        Update the time for the specified vertex based on global times
        """
        for e in self.v2e[v]:
            ix = self.elix(e, v)
            self.els[e].ts[ix] = self.ts[v]

    def gtdotns(self, v):
        return [np.round(self.els[e].gtdotn(self.elix(e, v)),3)
                for e in self.v2e[v]]

    def allgtdotns(self):
        return [self.gtdotns(v) for v in range(4)]

    def isoptimal(self, vtx):
        """
        a vertex is optimal if the
        gtdotn corresponding to its minimal kbar is positive.
        """
        gtdotns = self.gtdotns(vtx)
        ix = np.argmin(self.kbars(vtx))
        # print("vtx", vtx)
        # print("gtdotns", gtdotns)
        # print("self.kbars(vtx)", self.kbars(vtx))
        # print("ix", ix)
        return gtdotns[ix] > 0

    def kbars(self, v):
        """
        Given a vertex, compute time increments which would set ||∇ τ|| = 1/cbar
        for each adjacent element.
        """
        result = []
        for e in self.v2e[v]:
            el = self.els[e]
            eix = self.elix(e, v)
            result.append(el.kbar(eix))
        return result

    def ktildes(self):
        """
        Compute the time increment for each vertex v such that
        max over elements incident on v of ||∇ τ|| = 1/cbar
        """
        result = []
        for v in range(4):
            kbars = self.kbars(v)
            result.append(min(kbars))
        return result

    def weight(self, vtx):
        """
        function used to weight the ktilde for a vertex in the selection process
        """
        wt = 0.0
        for e in self.v2e[vtx]:
            vix = self.elix(e, vtx)
            term = self.els[e].Msums[vix]
            # print("vtx", vtx, "term", term)
            wt += term
        return wt

    def select_vertex(self, init2=0, init1=1, verbose=False):
        """
        select either v0 or v2, if one or the other is optimal.
        if both are optimal, select v0 for now.
        If neither is optimal, choose an optimal vertex in {v1, v3} which is
        minimal in time
        """
        if init2==0 and max(self.ts) == 0:
            if verbose:
                print("times 0")
            return 0
        elif init2==2 and max(self.ts) == 0:
            if verbose:
                print("times 0")
            return 2
        if self.isoptimal(0):
            if verbose:
                print("0 is optimal")
            return 0
        elif self.isoptimal(2):
            if verbose:
                print("2 is optimal")
            return 2
        if self.isoptimal(1) and self.isoptimal(3):
            if self.ts[1] == self.ts[3]:
                if verbose:
                    print("1 and 3 both optimal and equal")
                return init1
            elif self.ts[1] < self.ts[3]:
                if verbose:
                    print("1 and 3 both optimal; t(1) < t(3)")
                return 1
            else:
                if verbose:
                    print("1 and 3 both optimal; t(3) < t(1)")
                return 3
        elif self.isoptimal(1):
            if verbose:
                print("1 is optimal")
            return 1
        elif self.isoptimal(3):
            if verbose:
                print("3 is optimal")
            return 3
        else:
            raise ValueError("can't find a vertex")


    def step(self, verbose=False, weight=False, choose=False, init2=0, init1=1):
        """
        Select the vertex with maximal ktilde
        Pitch a tent at that vertex
        """
        ktildes = self.ktildes()
        if weight:
            prods = [k*self.weight(i) for i, k in enumerate(ktildes)]
            print("prods", prods)
            imax = np.argmax(prods)
        elif choose:
            if verbose:
                print("selecting a vertex")
            imax = self.select_vertex(init2=init2, init1=init1, verbose=verbose)
        else:
            imax = np.argmax(ktildes)
        ktmax = ktildes[imax]
        if ktmax < self.minstep:
            self.minstep = ktmax
        self.ts[imax] += ktmax
        self.updatets(imax)
        if verbose:
            # print("ktildes", ktildes)
            print("updated vertex", imax)
            # print("times", self.ts)
            print("minstep", self.minstep)
            # preview ktildes for next step
            # ktildes = self.ktildes()
            # print("next step ktildes", ktildes)
        return ktmax

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


def gen(verbose=False, weight=False, choose=False):
    ctr = 0
    restart = False
    inits = [[0,1],[2,1],[0,3],[2,3]]
    while ctr < 1000:
        retry_steps = []
        retries = 0
        while retries < 4:
            if restart:
                # restart the slab with the same vertices
                s = slab2d(vpts=s.vpts)
                if verbose:
                    print("restarting with same vertices")
                    print("times", s.ts, "vpts", s.vpts)
                restart = False
            else:
                # generate a new random slab
                if verbose:
                    print()
                    print()
                    print("Generating a new mesh")
                s = slab2d()
            for i in range(100):
                if verbose:
                    print("step", i)
                try:
                    ktmax = s.step(verbose=verbose, weight=weight, choose=choose,
                                   init2=inits[retries][0],
                                   init1=inits[retries][1])
                except(ValueError):
                    retry_steps.append(i)
                    retries += 1
                    if verbose:
                        print()
                        print("Can't find a vertex, retries = ", retries)
                        print()
                    restart = True
                    break # break out of for loop
            if not restart:
                break # break out of retry loop
        if retries >= 4:  # Houston: there is a problem
            print("counter", ctr)
            print("retry_steps", retry_steps)
            return s, retry_steps
        ctr += 1  # try another trig
    print("counter", ctr)

def multigen():
    while True:
        s, retry_steps = gen(choose=True)
        if sum(retry_steps) < 21:
            break
    s.plot()
    return s, retry_steps


# array([[0.2861798 , 0.59793209],
#        [0.63771079, 0.00240717],
#        [0.45567868, 0.95506259],
#        [0.40467456, 0.99532401]])
# 4 4 6 6 (not too bad)

# array([[0.2861798 , 0.59793209],
#        [0.63771079, 0.00240717],
#        [0.45567868, 0.95506259],
#        [0.40467456, 0.99532401]])
# 4 4 7 6 (not too bad)

# array([[0.58228029, 0.18701104],
#        [0.99382862, 0.81330621],
#        [0.26723904, 0.27925846],
#        [0.22784489, 0.05875047]])
# 4 5 7 7 (nice looking)

# array([[0.02838083, 0.78813412],
#        [0.60655684, 0.23317332],
#        [0.08799146, 0.99574778],
#        [0.03615092, 0.99614195]])
# 5 6 5 6

if __name__ == "__main__":

    generate = False
    example = 4

    if generate:
        s = gen()
    else:
        if example == 1:
            vpts = np.array([[2, 0], [6.0, 0], [3, np.sqrt(3)], [0, 0]])
        elif example == 2:
            vpts = np.array([[0, 0], [5, 1.1], [0, 1.], [-3., -4.]])
        elif example == 3:
            vpts = np.array([[1.2, 0], [1.9, .5], [.7, .5], [0, 0]])
        elif example == 4:
            vpts = np.array([[0.55407862, 0.24615304],
                             [0.72966765, 0.503949],
                             [0.53244779, 0.27788538],
                             [0.53841842, 0.05414372]])
        s = slab2d(vpts=vpts)
