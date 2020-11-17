import matplotlib.pyplot as plt
import ngsolve as ng
import numpy as np

from .tents import (TentPitchedSlab1, TentPitchedSlab2, TentPitchedSlab3)


class TentSlab(object):
    def __init__(self, mesh, method=None, heapsize=None):
        if method is None:
            method = "edge"
        self.mesh = mesh
        self.dim = mesh.dim
        if self.dim == 1:
            if heapsize is None:
                self.slab = TentPitchedSlab1(mesh, method)
            else:
                self.slab = TentPitchedSlab1(mesh, method,
                                             heapsize=heapsize)
        elif self.dim == 2:
            if heapsize is None:
                self.slab = TentPitchedSlab2(mesh, method)
            else:
                self.slab = TentPitchedSlab2(mesh, method,
                                             heapsize=heapsize)
        elif self.dim == 3:
            if heapsize is None:
                self.slab = TentPitchedSlab3(mesh, method)
            else:
                self.slab = TentPitchedSlab3(mesh, method,
                                             heapsize=heapsize)
        else:
            raise NotImplementedError("mesh dimension not supported")

    def SetWavespeed(self, c):
        self.slab.SetWavespeed(c)

    def PitchTents(self, dt, local_ct=None, global_ct=None):
        print("****************\n"
              "IMPORTANT NOTICE\n"
              "****************\n"
              "In case the slab can not be pitched"
              " execution will NOT be aborted.\n"
              "The return value of this function must ALWAYS be checked.\n"
              "Pitching slab...")
        if local_ct is None:
            local_ct = False
        if global_ct is None:
            success = self.slab.PitchTents(dt, local_ct)
        else:
            success = self.slab.PitchTents(dt, local_ct, global_ct)
        if success is True:
            print("The slab was successfully pitched!")
        else:
            print("The slab could not be pitched.\n"
                  "If desired, it can still be printed"
                  " for debugging purposes")
        return success

    def GetNTents(self):
        return self.slab.GetNTents()

    def GetNLayers(self):
        return self.slab.GetNLayers()

    def GetSlabHeight(self):
        return self.slab.GetSlabHeight()

    def MaxSlope(self):
        return self.slab.MaxSlope()

    def DrawPitchedTentsVTK(self):
        if self.dim != 2:
            raise NotImplementedError("Only supported for 2D spatial mesh")
        return self.slab.DrawPitchedTentsVTK()

    def DrawPitchedTentsGL(self):
        if self.dim == 1:
            raise NotImplementedError("1D spatial mesh not supported")
        elif self.dim == 2:
            return self.slab.DrawPitchedTentsGL()
        else:
            # times are not used in 3D case
            data, ntents, nlevels = self.slab.DrawPitchedTentsGL()
            return data, None, ntents, nlevels

    def DrawPitchedTentsPlt(self):
        """
        Draw 1D tents using MatPlotlib
        """
        if self.dim != 1:
            raise NotImplementedError("Only supported for 1D spatial mesh")
        else:
            tents = self.slab.DrawPitchedTentsPlt()
            pnts = []
            colors = list('bgrcmyk')
            for pnt in self.mesh.ngmesh.Points():
                pnts.append(pnt.p[0])
            for i, tent in enumerate(tents):
                layer = tent[0][3] % 7
                if len(tent) == 3:
                    xvals = [pnts[tent[1][0]],
                             pnts[tent[0][0]], pnts[tent[2][0]]]
                    tvals = [tent[1][1], tent[0][1], tent[2][1]]
                    xpos = pnts[tent[0][0]]
                    tpos = 1/3*(tent[1][1]+tent[0][1]+tent[2][1])
                else:
                    xvals = [pnts[tent[1][0]], pnts[tent[0][0]]]
                    tvals = [tent[1][1], tent[0][1]]
                    xpos = 0.5*(pnts[tent[1][0]]+pnts[tent[0][0]])
                    tpos = 0.5*(tent[1][1]+tent[0][1])
                plt.plot(xvals, tvals, color=colors[layer])
                plt.text(xpos, tpos, str(i), horizontalalignment='center',
                         verticalalignment='center')
            plt.ylim([0, self.dt*1.1])
            plt.show()

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

        mvs = [mesh[ng.NodeId(ng.VERTEX, v)] for v in [vtx]+nbv]
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

    def GetTent(self, nr):
        return self.slab.GetTent(nr)
