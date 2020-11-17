from .tents import (Burgers1, Burgers2)


class Burgers(object):
    def __init__(self, tentslab, order, flags={}):
        """
        INPUTS:

        tentslab : TentSlab
        order    : int : order of the L2 space
        flags    : optional dict of NGSolve flags
        """
        self.slab = tentslab.slab
        self.order = order
        self.flags = flags
        self.dim = tentslab.dim
        if self.dim == 1:
            self.conslaw = Burgers1(self.slab, self.order, self.flags)
        elif self.dim == 2:
            self.conslaw = Burgers2(self.slab, self.order, self.flags)
        else:
            raise NotImplementedError("Burgers not implemented for 3D")
        self.sol = self.conslaw.sol
        self.res = self.conslaw.res
        self.nu = self.conslaw.nu
        self.space = self.conslaw.space

    def SetInitial(self, cf):
        self.conslaw.SetInitial(cf)

    def PropagatePicard(self, vecu, steps):
        self.conslaw.PropagatePicard(vecu, steps)

    def Tau(self):
        return self.conslaw.Tau()
