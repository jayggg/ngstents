from ._pyconslaw import ConservationLaw
from ._pyconslaw import ApplyDG

class Burgers(ConservationLaw):
    def __init__(self, gfu, tentslab):
        """
        INPUTS:
        
        gfu      : GridFunction
        tentslab : TentSlab
        """
        ConservationLaw.__init__(self, gfu, tentslab, "burgers")

class Euler(ConservationLaw):
    def __init__(self, gfu, tentslab):
        """
        INPUTS:
        
        gfu      : GridFunction
        tentslab : TentSlab
        """
        ConservationLaw.__init__(self, gfu, tentslab, "euler")

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
