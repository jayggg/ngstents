from .conslaw import ConservationLaw

def Burgers(tentslab, order):
    """
    INPUTS:
    
    tentslab : TentSlab
    order    : int : order of the L2 space
    """
    return ConservationLaw(tentslab,"burgers",order)


def Wave(tentslab, order):
    """
    INPUTS:
    
    tentslab : TentSlab
    order    : int : order of the L2 space
    """
    return ConservationLaw(tentslab, "wave", order)

# clean up
del conslaw
