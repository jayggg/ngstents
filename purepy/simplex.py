import numpy as np
import numpy.linalg as la


class simplex(object):
    """
    A simplex on which we build a tent slab
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
        self.ts = None
        Ms, Msums, Mdot1s, Ainvs = [], [], [], []
        for i in range(self.d+1):
            At = np.roll(pts - pts[i, :], -(i+1), 0)[:self.d, :]
            Ainvs.append(la.inv(At.T))
            M = la.inv(At@At.T)
            Ms.append(M)
            Mdot1s.append(M.sum(axis=0))
            Msums.append(M.sum())
        self.Ms, self.Msums, self.Mdot1s, self.Ainvs = Ms, Msums, Mdot1s, Ainvs

    def check(self, ts, k):
        """
        Check if the condition ∇ τ ⋅ nₖ > 0
        holds for the vertex k

        INPUTS

        ts: numpy array with 4 elements representing vertex times
        k: vertex index to check
        """
        return self.gradtau(ts).dot(self.normals()[k]) >= 0

    def normals(self):
        """
        Compute the normal vectors to each facet indexed by opposite vertex
        """
        if self.nmls is None:
            # normal vector of facet opposite v0 for unit tet
            nml0 = np.ones(self.d)
            nmls = []
            for i in range(self.d + 1):
                nmls.append(nml0.dot(self.Ainvs[i]))
            self.nmls = np.array(nmls)
        return self.nmls

    def gradtau(self, ts):
        """
        Compute ∇ τ for the given vertex times

        INPUTS

        ts: numpy array with d+1 elements representing vertex times
        """
        dts = self.dts(ts, 0)
        return dts.dot(self.Ainvs[0])

    def gradnormsq(self, ts):
        """
        Compute ||∇ τ||²

        INPUTS

        ts: numpy array with d+1 elements representing vertex times
        """
        dts = self.dts(ts, 0)
        return dts.dot(self.Ms[0]).dot(dts)

    def dts(self, ts, k):
        """
        Compute dtₖ, the vector of time differences relative to vₖ

        INPUTS

        ts: numpy array with d+1 elements representing vertex times
        k: index of reference vertex
        """
        return np.roll(ts - ts[k], -(k+1))[:self.d]

    def ktilde(self, ts, k):
        """
        Compute k̃ₖ for vertex k

        INPUTS

        ts: numpy array with d+1 elements representing vertex times
        """
        M, Msum, Mdot1 = self.Ms[k], self.Msums[k], self.Mdot1s[k]
        c = self.cbar
        dt = self.dts(ts, k)
        m1dt = Mdot1.dot(dt)
        disc = m1dt**2 - Msum * (dt.dot(M).dot(dt) - 1/c)
        if disc < 0:
            disc = 0.
        return (m1dt + np.sqrt(disc))/Msum

    def gradangles(self, ts):
        """
        Compute the spherical angles φ and θ for ∇ τ in 3D
        or the angle θ in 2D

        INPUTS
        ts: numpy array with d+1 elements representing vertex times
        """
        if self.d not in [2, 3]:
            # angles only supported for 2D and 3D
            return None
        grad = self.gradtau(ts)
        theta = np.arctan2(grad[1], grad[0])*180/np.pi
        if self.d == 3:
            ngrad = la.norm(grad)
            if ngrad == 0:
                return (0, 0)
            phi = np.arccos(grad[2]/ngrad)*180/np.pi
            return (phi, theta)
        else:
            return theta

    def nmlangles(self):
        """
        Compute spherical angles φ and θ of the opposite normal vectors
        in 3D or the angle θ in 3D
        """
        if self.d not in [2, 3]:
            # angles only supported for 2D and 3D
            return np.zeros((self.d+1, 2))
        nmls = self.normals()
        theta = np.arctan2(nmls[:, 1], nmls[:, 0])*180/np.pi
        if self.d == 3:
            nnmls = ((nmls**2).sum(axis=1))**(1./2)
            phi = np.arccos(nmls[:, 2]/nnmls)*180/np.pi
            return np.array(list(zip(phi, theta)))
        else:
            return theta

    def dists(self):
        """
        Compute distance from each vertex to its opposite face
        This should equal sqrt(1/self.Msums)
        """
        pts = self.pts
        V = np.roll(pts, -1, 0)-pts
        vol = np.linalg.det(V[:self.d, :])
        areas = (((np.cross(V, np.roll(V, -1, 0)))**2).sum(1))**(.5)
        return vol/np.roll(areas, -1)

    def iterate(self, t0=None):
        """
        Pitch tents, printing information and checking gradient condition

        INPUTS
        t0: optional numpy array with 4 elements representing vertex times
        """
        def printangles():
            if self.d == 3:
                print("gradtau: phi: {:.1f} theta: {:.1f}".format(*ga))
                print("oppnml:  phi: {:.1f} theta: {:.1f}".format(*na))
            elif self.d == 2:
                print("gradtau: theta: {:.1f}".format(ga))
                print("oppnml:  theta: {:.1f}".format(na))

        if t0 is None:
            t0 = np.zeros(self.d+1)
        self.ts = t0
        while True:
            ktildes = [self.ktilde(self.ts, i) for i in range(self.d+1)]
            mx = max(ktildes)
            imx = np.argmax(ktildes)
            ga = self.gradangles(self.ts)
            na = self.nmlangles()[imx, :]
            print("ktildes", ktildes)
            print("mx", mx)
            print("imx", imx)
            printangles()
            if self.check(self.ts, imx):
                print("∇ τ ⋅ n > 0 for vertex with max k̃")
            else:
                raise ValueError("Expected gradient condition to hold")
            self.ts[imx] += mx
            print("self.ts", self.ts)
            resp = input("press enter or q ")
            if resp.lower() == 'q':
                break
            print()


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
    pts = np.zeros((d+1, d))
    pts[1:, :] = np.eye(d)
    return simplex(pts, cbar)
