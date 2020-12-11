from .conslaw import ConservationLaw

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

# clean up
del conslaw
