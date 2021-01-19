from .conslaw import ConservationLaw
from .conslaw import ApplyDG

class Burgers(ConservationLaw):
    def __init__(self, tentslab, order):
        """
        INPUTS:
        
        tentslab : TentSlab
        order    : int : order of the L2 space
        """
        ConservationLaw.__init__(self, tentslab, "burgers", order)

class Wave(ConservationLaw):
    def __init__(self, tentslab, order):
        """
        INPUTS:
        
        tentslab : TentSlab
        order    : int : order of the L2 space
        """
        ConservationLaw.__init__(self, tentslab, "wave", order)

class Advection(ConservationLaw):
    def __init__(self, tentslab, order):
        """
        INPUTS:
        
        tentslab : TentSlab
        order    : int : order of the L2 space
        """
        ConservationLaw.__init__(self, tentslab, "advection", order)

class Maxwell(ConservationLaw):
    def __init__(self, tentslab, order):
        """
        INPUTS:
        
        tentslab : TentSlab
        order    : int : order of the L2 space
        """
        ConservationLaw.__init__(self, tentslab, "maxwell", order)
        
# clean up
del conslaw
