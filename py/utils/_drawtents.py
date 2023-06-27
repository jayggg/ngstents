import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def DrawPitchedTentsPlt(self, showtentnums=False):
    """
    Draw 1D tents using Matplotlib
    """
    if self.mesh.dim != 1:
        raise NotImplementedError("Only supported for 1D spatial mesh")
    else:
        tents = self.TentData1D()
        pnts = []
        fig, ax = plt.subplots()

        colors = list(mcolors.BASE_COLORS)[:-1]  # remove white
        colors += list(mcolors.TABLEAU_COLORS)[:-1]  # omit tab:cyan close to c
        for pnt in self.mesh.ngmesh.Points():
            pnts.append(pnt.p[0])
        for i, tent in enumerate(tents):
            layer = tent[0][3] % len(colors)
            if len(tent) == 3:
                xvals = [pnts[tent[1][0]], pnts[tent[0][0]], pnts[tent[2][0]]]
                tvals = [tent[1][1], tent[0][1], tent[2][1]]
                xpos = pnts[tent[0][0]]
                tpos = 1 / 3 * (tent[1][1] + tent[0][1] + tent[2][1])
            else:
                xvals = [pnts[tent[1][0]], pnts[tent[0][0]]]
                tvals = [tent[1][1], tent[0][1]]
                xpos = 0.5 * (pnts[tent[1][0]] + pnts[tent[0][0]])
                tpos = 0.5 * (tent[1][1] + tent[0][1])
            ax.plot(xvals, tvals, color=colors[layer])
            if showtentnums:
                ax.text(xpos,
                        tpos,
                        str(i),
                        horizontalalignment='center',
                        verticalalignment='center')

        ax.set_ylabel('time ($t$)')
        ax.set_xlabel('space ($x$)')
        ax.margins(x=0, y=0)
        ax.axis('scaled')
        plt.show()
