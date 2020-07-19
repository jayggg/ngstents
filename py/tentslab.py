import matplotlib.pyplot as plt

from tents import (TentPitchedSlab1, TentPitchedSlab2, TentPitchedSlab3)


class TentSlab(object):
    def __init__(self, mesh, dt, c, heapsize=None):
        self.mesh = mesh
        self.dt = dt
        self.c = c
        self.dim = mesh.dim
        if self.dim == 1:
            if heapsize is None:
                self.slab = TentPitchedSlab1(mesh, dt, c)
            else:
                self.slab = TentPitchedSlab1(mesh, dt, c, heapsize=heapsize)
        elif self.dim == 2:
            if heapsize is None:
                self.slab = TentPitchedSlab2(mesh, dt, c)
            else:
                self.slab = TentPitchedSlab2(mesh, dt, c, heapsize=heapsize)
        elif self.dim == 3:
            if heapsize is None:
                self.slab = TentPitchedSlab3(mesh, dt, c)
            else:
                self.slab = TentPitchedSlab3(mesh, dt, c, heapsize=heapsize)
        else:
            raise NotImplementedError("mesh dimension not supported")

    def GetNTents(self):
        return self.slab.GetNTents()

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


    def GetTent(self, nr):
        return self.slab.GetTent(nr)
