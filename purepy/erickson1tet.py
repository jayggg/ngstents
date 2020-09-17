import numpy as np
import numpy.linalg as la


class Tet1mesh(object):
    """
    A tent slab consisting of one tetrahedron
    """

    def __init__(self, pts, cbar=1.0, eps1=.1, eps2=.6):
        """
        INPUTS

        pts: 4x3 numpy array, each row a vertex of a tetrahedron
        cbar: float wavespeed
        eps1: progress constraint constant
        eps2: progress constraint constant based on 2D cone constraint
        """
        self.pts = pts
        self.cbar = cbar
        self.eps1 = eps1
        self.eps2 = eps2
        self.ts = np.zeros(4)  # vertex times

    def oppfacet(self, vix):
        """
        return the vertices of the opposite facet to a vertex
        """
        i, j, k = (vix+1) % 4, (vix+2) % 4, (vix+3) % 4
        return [i, j, k]

    def trigswithvtx(self, vix):
        """
        return a list of pairs, with the first element of each pair a list
        representing a triangle containing the given vertex, placing
        the given vertex first in the list, and with the second element
        of the pair the remaining vertex (not in the triangle)
        """
        i, j, k = (vix+1) % 4, (vix+2) % 4, (vix+3) % 4
        return [([vix, j, k], i), ([vix, k, i], j), ([vix, i, j], k)]

    def timeorderedpts(self, ixs):
        """
        return vertices ordered with ascending time values
        """
        _, tq, tr = self.ts[ixs]
        if tr < tq:
            ixs[1], ixs[2] = ixs[2], ixs[1]
        return self.pts[ixs, :], self.ts[ixs]

    def stpts(self, ixs=None):
        """
        return spacetime points (optionally for specified vertices)
        """
        if ixs is None:
            pts = self.pts
            ts = self.ts
        else:
            pts = self.pts[ixs]
            ts = self.ts[ixs]
        return np.hstack((pts, ts[:, np.newaxis]))

    def progconstr(self, ixs, ix):
        """
        return the progress constraint for a given trig and point outside it
        """
        pts, ts = self.timeorderedpts(ixs)
        tp, tq, tr = ts
        p, q, r = pts
        po = self.pts[ix, :]
        pq, rq, trq = p-q, r-q, tr-tq
        nrq = la.norm(rq)
        wp = la.norm(np.cross(pq, rq)) / nrq
        sigmaF = calc_sigmaF(pts, po)
        # set conetp so that ∥∇_H τ∥ = (1-ε₂) σ_F / c̅
        disc = (nrq*sigmaF*(1-self.eps2)/self.cbar)**2 - trq**2
        if disc < 0:
            print("negative disc", disc, "setting to zero")
            disc = 0
        conetp = tq + trq*pq.dot(rq)/nrq**2 + wp*np.sqrt(disc)/nrq
        print("2D cone progress constraint: {:.3f}".format(conetp-tp))
        prog_constr = tr + (1-self.eps1)*wp/self.cbar
        print("2D progress constraint: {:.3f}".format(prog_constr-tp))
        newtp = min(conetp, prog_constr)
        return newtp - tp

    def coneconstr(self, vix):
        """
        return the global cone constraint for a given vertex
        """
        ixs = self.oppfacet(vix)
        pts = self.pts[ixs]
        stpts = self.stpts(ixs)
        p = self.pts[vix]
        tp = self.ts[vix]
        pH = calc_pH(pts, p)
        pHt = calc_pHt(stpts, pH)
        sigmaF = calc_sigmaF(pts, p)
        term = self.eps2 + (1-self.eps2)*np.sqrt(1-sigmaF**2)
        newtp = pHt + la.norm(p-pH)*term/self.cbar
        print("3D cone constraint: {:.3f}".format(newtp-tp))
        return newtp - tp

    def ktilde(self, vix):
        """
        Compute a time increment k̃ for a minimal vertex
        """
        dts = []
        for ixs, ix in self.trigswithvtx(vix):
            dts.append(self.progconstr(ixs, ix))
        dts.append(self.coneconstr(vix))
        return max(0, min(dts))

    def iterate(self):
        """
        Pitch tents on the tetraheron, using the Erickson higher dimensional
        algorithm, printing information and checking gradient condition
        """
        msg = "constraints prevented lifting vertex to ε₂∥p-pH∥ = {:.3f}!"
        print("pitching tents on a 1-tet mesh with eps1 = {}, eps2 = {}\n"
              .format(self.eps1, self.eps2))
        while True:
            vix = np.argmin(self.ts)
            print("minimal vertex:", vix)
            ktilde = self.ktilde(vix)
            print("ktilde: {:.3f}".format(ktilde))
            self.ts[vix] += ktilde
            print("self.ts", self.ts)
            self.check()
            newmin = min(self.ts)
            ixs = self.oppfacet(vix)
            pts = self.pts[self.oppfacet(vix)]
            p = self.pts[vix]
            pH = calc_pH(pts, p)
            mintp = newmin + self.eps2*la.norm(p-pH)
            if self.ts[vix] < mintp:
                print(msg.format(mintp))
            if ktilde < 1.0e-15:
                print("ktilde zero - degenerate tent. Stopping.")
                break
            resp = input("press enter (or q to quit)")
            if resp.lower() == 'q':
                break
            print()

    def check(self):
        """
        verify that all facets meet satisfy the progress constraint
        and that the tet satisfies the global cone constraint
        """
        tol = 1.0e-14
        for i in range(4):
            ixs = self.oppfacet(i)
            gtnorm = calc_gradHtaunorm(self.stpts(ixs))
            sigmaF = calc_sigmaF(self.pts[ixs], self.pts[i])
            condition = gtnorm <= (1-self.eps2)*sigmaF/self.cbar + tol
            if not condition:
                print("gtnorm", gtnorm)
                print("(1-eps2)*sigmaF/cbar", (1-self.eps2)*sigmaF/self.cbar)
            assert condition
        normgrad = calc_normgrad3D(self.stpts())
        assert normgrad <= 1/self.cbar

    def plot(self):
        """
        Scatter plot of the vertices of the tetrahedron
        """
        import matplotlib.pyplot as plt
        pts = self.pts
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2])
        plt.show()


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

    assert sigmaF <= 1 + 1.0e-15
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

    p: a 3D point, not necessarily in the plane containing the triangle F
    v0, v1, v2: vertices of triangle F
    w0, w1, w2: vectors def by w0 = v1-v0, w1 = v2-v1, w2 = v0-v2
    w0sq, w1sq, w2sq: squared norms of vectors
    """
    p = np.array(p)
    v0, v1, v2 = A

    w10, w21, w02 = np.roll(A, -1, 0) - A
    w20, w01, w12 = -w02, -w10, -w21

    w01sq, w12sq, w20sq = w01.dot(w01), w12.dot(w12), w20.dot(w20)

    # project p onto plane of triangle with verts in A
    pH = calc_pH(A, p)

    p0, c0 = closest_pt(pH, v0, v2, w10, w20, w20sq)
    p1, c1 = closest_pt(pH, v1, v0, w21, w01, w01sq)
    p2, c2 = closest_pt(pH, v2, v1, w02, w12, w12sq)

    if p0 is None and p1 is None and p2 is None:
        return pH   # inside triangle
    if p0 is not None and 0 <= c0 <= 1:
        return p0   # on edge w20
    if p1 is not None and 0 <= c1 <= 1:
        return p1   # on edge w01
    if p2 is not None and 0 <= c2 <= 1:
        return p2   # on edge w12
    if p0 is not None:
        return p0   # v0 if c0 < 0 else v2
    if p1 is not None:
        return p1   # v1 if c1 < 0 else v0
    return p2       # v2 if c2 < 0 else v1


def closest_pt(p, v0, v2, w10, w20, w20sq):
    """
    INPUTS

    p, v0, v1: 3D points in same plane
    w10 = v1-v0
    w20 = v2-v0
    w20sq = norm squared of w20
    """
    # check if p and v1 on same side of w20
    d = np.cross(p-v0, w20).dot(np.cross(w10, w20))
    if d <= 0:
        # opposite sides so check projection
        c = (p-v0).dot(w20)/w20sq
        if c < 0:
            return v0, c
        elif c > 1:
            return v2, c
        return v0 + c*w20, c  # projection
    else:
        return None, None


def calc_normgrad3D(pts):
    """
    Given a list of four spacetime points, compute the normal
    to the hyperplane they lie in and its gradient
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
    Return a matrix of point differences and minor matices
    """
    V = A[1:, :] - A[:-1, :]
    return V[:, [1, 2, 3]], V[:, [0, 2, 3]], \
        V[:, [0, 1, 3]], V[:, [0, 1, 2]], V


def rand_simplex(eps1=0.1, eps2=0.5):
    """
    Generate a random d-simplex with vertices ordered so that it is
    homeomorphic to the unit d-simplex
    """
    pts = np.random.rand(4, 3)
    return Tet1mesh(pts, cbar=1, eps1=eps1, eps2=eps2)


def atet(eps1=0.1, eps2=0.3):
    """
    A tetrahedron containing an obtuse angle.
    Tent pitching fails for this element if eps2 >= .43
    """
    pts = np.array([[1, .333, 0], [.866, .5, .127], [1, .373, .127],
                    [1, .54, .127]])
    return Tet1mesh(pts, cbar=1, eps1=eps1, eps2=eps2)


np.set_printoptions(precision=3)
s = atet(eps1=.1, eps2=.37)
s.iterate()
