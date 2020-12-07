from .conslaw import ConservationLaw

def Burgers(tentslab, order):
    """
    INPUTS:
    
    tentslab : TentSlab
    order    : int : order of the L2 space
    """
    return ConservationLaw(tentslab,"burgers",order)

# clean up
del conslaw
