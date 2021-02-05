from .conslaw import ConservationLaw
from .conslaw import ApplyDG

# class Burgers(ConservationLaw):
#     def __init__(self, tentslab, order):
#         """
#         INPUTS:
#         
#         tentslab : TentSlab
#         order    : int : order of the L2 space
#         """
#         ConservationLaw.__init__(self, tentslab, "burgers", order)
# 
class Wave(ConservationLaw):
    def __init__(self, gfu, tentslab):
        """
        INPUTS:
        
        gfu      : GridFunction
        tentslab : TentSlab
        """
        ConservationLaw.__init__(self, gfu, tentslab, "wave")

class Advection(ConservationLaw):
    def __init__(self, gfu, tentslab):
        """
        INPUTS:
        
        gfu      : GridFunction
        tentslab : TentSlab
        """
        ConservationLaw.__init__(self, gfu, tentslab, "advection")

class Maxwell(ConservationLaw):
    def __init__(self, gfu, tentslab):
        """
        INPUTS:
        
        gfu      : GridFunction
        tentslab : TentSlab
        """
        ConservationLaw.__init__(self, gfu, tentslab, "maxwell")
        
# clean up
del conslaw
