import matplotlib.pyplot
import ngsolve
import numpy

def DrawPitchedTentsPlt(self):
    """
    Draw 1D tents using MatPlotlib
    """
    if self.mesh.dim != 1:
        raise NotImplementedError("Only supported for 1D spatial mesh")
    else:
        # tents = self.DrawPitchedTents()
        tents = self._TentData1D()
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
            matplotlib.pyplot.plot(xvals, tvals, color=colors[layer])
            matplotlib.pyplot.text(xpos, tpos, str(i), horizontalalignment='center',
                     verticalalignment='center')
        matplotlib.pyplot.ylim([0, self.GetSlabHeight()*1.1])
        matplotlib.pyplot.show()

def Draw3DTentPlt(self, tentnr):
    """
    Draw a single 3D tent using Matplotlib. Colors are used to represent
    the times of the neighbor vertices.  The tent pole height is represented
    by the size of the central vertex.
    """
    if self.mesh.dim != 3:
        raise NotImplementedError("Only supported for 3D spatial mesh")
    mesh = self.mesh
    tent = self.GetTent(tentnr)
    vtx = tent.vertex
    nbv = list(tent.nbv)
    tt = tent.ttop
    tb = tent.tbot
    nbtime = list(tent.nbtime)

    mvs = [mesh[ngsolve.NodeId(ngsolve.VERTEX, v)] for v in [vtx]+nbv]
    mels = [mesh[ngsolve.ElementId(ngsolve.VOL, e)] for e in tent.els]
    pts = numpy.array([v.point for v in mvs])
    facetvs = [[mvs.index(v) for v in el.vertices if v != mvs[0]]
               for el in mels]
    fig = matplotlib.pyplot.figure()
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
          .format(tt, tb, tt-tb))
    return ax

# TentSlab.DrawPitchedTentsPlt = DrawPitchedTentsPlt
# TentSlab.Draw3DTentPlt = Draw3DTentPlt
# del DrawPitchedTentsPlt
# del Draw3DTentPlt
