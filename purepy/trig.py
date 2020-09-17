import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la


class TriangleSlab(object):
    def __init__(self, pts, t=None, cbar=None, slabheight=None):
        """
        INPUTS

        pts: a 3x2 numpy array containing the three 2D vertices
        cbar: wavespeed (a positive float)
        """
        # 2D vertices
        self.pts = pts
        # the wavespeed
        self.cbar = 1.0 if cbar is None else cbar
        # the height (in time) of the desired slab
        self.slabheight = 1.0 if slabheight is None else slabheight
        # rolling differences of vertices
        self.w = np.roll(pts, -1, 0) - pts
        # norms squared of rows of v
        self.wsq = (self.w**2).sum(1)
        # dot products of successive rows of v, i.e. v0⋅v1, v1⋅v2, v2⋅v0
        self.wdot = (self.w*np.roll(self.w, -1, 0)).sum(1)
        # the square of spatial determinant/wavespeed
        self.csq = (la.det(self.w[:2, :2])/cbar)**2
        # array of vertex times
        self.t = np.zeros(3) if t is None else t
        # array of vertex time diffs: t(v1)-t(v0), t(v2)-t(v1), t(v0)-t(v2)
        self.update_dt()
        # previous tent vertex index for testing
        self.previx = None

    def A(self):
        """
        Get a 3x3 matrix with spacetime points
        """
        return np.hstack((self.pts, self.t[:, np.newaxis]))

    def V(self):
        """
        Get a 3x3 matrix with rows the cyclic differences of the
        rows of A
        """
        A = self.A()
        return np.roll(A, -1, 0) - A

    def gradnorm(self):
        """
        Compute the norm of the gradient of the time function tau
        """
        V = self.V()
        nml = np.cross(V[0, :], V[1, :])
        return ((nml[:2]**2).sum()/nml[2]**2)**(1/2.)

    def gradopt(self):
        """
        Determine which vertices will have maximal ktilde based on the
        direction of the gradient relative to the opposite edge vectors.
        If the y component of the difference is negative the vertex
        should have maximal ktilde.
        """
        V = self.V()
        nml = np.cross(V[0, :], V[1, :])
        grad = -nml[:2]/nml[2]
        return np.cross(grad, np.roll(self.w, -1, 0)) >= 0

    def checkdt(self, i, newdt):
        """
        Pass me the index of the vertex which has the largest computed
        ktilde value.  I will check if the gradient has optimal direction
        for this vertex.
        """
        gopt = self.gradopt()[i]
        print("gradient direction is{} optimal for v{}"
              .format("" if gopt else " NOT", i))

    def getdelt(self, i):
        """
        Get a dt value for vertex_i which if added to tau(vertex_i)
        will result in a tent with gradient norm equal 1/wavespeed.
        Return 0 if no such dt value exists for the vertex
        """
        j = (i+1) % 3
        dti, dtj = self.dt[i], self.dt[j]
        vsq, vdot, csq = self.wsq, self.wdot, self.csq
        disc = (vdot[i]*dtj)**2 - vsq[j]*(vsq[i]*dtj*dtj - csq)
        print("disc", disc)
        term = np.sqrt(disc) if disc >= 0 else 0
        delti = min(self.slabheight - self.t[i],
                    dti - (vdot[i]*dtj - term)/vsq[j])
        return delti

    def maxdeltix(self, verbose=False):
        """
        Get the index and delt of the vertex with maximal computed delt value
        """
        delts = []
        for i in range(3):
            delts.append(self.getdelt(i))
        ix = np.argmax(delts)
        if verbose:
            print("delts", delts)
        return ix, delts

    def update_dt(self):
        """
        Update the array of vertex time differences
        """
        self.dt = np.roll(self.t, -1) - self.t

    def trypitch(self, verbose=False):
        """
        try to pitch a tent with maximal ktilde
        """
        ix, delts = self.maxdeltix(verbose=True)
        delt = delts[ix]
        if self.t[ix] + delt < self.slabheight:
            oix = (ix+1) % 3
            if oix != self.previx:
                err = self.checkdt(oix, delts[oix])
            oix = (ix+2) % 3
            if oix != self.previx:
                err = self.checkdt(oix, delts[oix])
            self.checkdt(ix, delt)
        self.previx = ix
        if delt > 0:
            self.t[ix] += delt
            self.update_dt()
            times = ", ".join(["{:.3f}".format(t) for t in self.t])
            print("Chose vertex {} with maximal delt {:.3f}; times are {}"
                  .format(ix, delt, times))
            print("New gradient norm: {:.3f}".format(self.gradnorm()))
        else:
            raise ValueError("No time increment for vertex with max ktilde")
        if verbose:
            print("A", self.A())
            print("V", self.V())
        err = abs(self.gradnorm() - 1/self.cbar)
        if err > 1.e-12 and self.t[ix] != self.slabheight:
            print("error in gradient {:.3e} exceeds tolerance".format(err))
        return False

    def maxedge(self):
        """
        return the length of the maximum edge of the triangle
        """
        return np.sqrt(max(self.wsq))

    def maxangles(self):
        """
        return the largest, then the next largest angles of the trig
        """
        cosines = -self.wdot/np.sqrt(self.wsq * np.roll(self.wsq, -1, 0))
        angles = sorted(np.arccos(cosines) * 180/np.pi, reverse=True)
        return angles[:2]

    def complete(self):
        """
        condition for complete tent pitched slab
        """
        return np.all(self.t >= self.slabheight)

# ----------------
# Helper functions
# ----------------


def plotpts(points):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(points[:, 0], points[:, 1], 'ro')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    ax.set_aspect('equal', adjustable='box')
    plt.show()


def plotdata(data, idx):
    points = data[idx][3]
    plotpts(points)


def gendata(cbar=1.0, slabheight=1.0):
    angles = []
    skews = []
    tents = []
    allpts = []
    for i in range(1000):
        pts = np.random.rand(3, 2)
        t = TriangleSlab(pts, cbar, slabheight)
        ntents = 0
        while not t.complete():
            # ix, delts = t.maxdeltix()
            err = t.trypitch()
            ntents += 1
            if ntents > 5000 or err:
                break
        bigangle, nextangle = t.maxangles()
        isoc = (180-bigangle)/2
        skew = abs(isoc-nextangle)/isoc
        print("vertices:", pts)
        print("times:", t.t)
        print("tents:", ntents)
        print("max edge length: {:.3f}".format(t.maxedge()))
        print("max angle: {:.3f} degrees".format(bigangle))
        print("skew: {:.3f}".format(skew))
        tents.append(ntents)
        angles.append(bigangle)
        skews.append(skew)
        allpts.append(pts)
        if err:
            break

    nangles = (np.array(angles)-60)/120
    plt.scatter(nangles, tents)
    plt.show()

    data = list(zip(np.arange(len(tents)), nangles, tents, allpts))
    data.sort(key=lambda d: d[1])


def checktrig(pts=None, t=None, cbar=1.0, slabheight=1.0):
    """
    Check the tent slab for a provided trig.  If not provided, a random
    trig will be used.
    """
    maxtents = 10  # 500
    if pts is None:
        pts = np.random.rand(3, 2)
    if t is None:
        t = np.zeros(3)
    w = np.roll(pts, -1, 0) - pts
    if np.cross(w[0, :], w[1, :]) < 0:
        print("forcing ccw ordering of vertices")
        pts[:2, :] = pts[1::-1, :]
    t = TriangleSlab(pts, t, cbar, slabheight)
    print("initial gradient norm", t.gradnorm())
    ntents = 0
    while not t.complete():
        t.trypitch()
        ntents += 1
        if ntents > maxtents:
            print("more than {} tents -- quitting".format(maxtents))
            break
    bigangle, nextangle = t.maxangles()
    isoc = (180-bigangle)/2
    skew = abs(isoc-nextangle)/isoc
    print("vertices:", pts)
    print("times:", t.t)
    print("tents:", ntents)
    print("max edge length: {:.3f}".format(t.maxedge()))
    print("max angle: {:.3f} degrees".format(bigangle))
    print("skew: {:.3f}".format(skew))
