from ._pyconslaw import ConservationLaw


class Burgers(ConservationLaw):
    def __init__(self, gfu, tentslab, **kwargs):
        """
        INPUTS:

        gridfunction : GridFunction
        tentslab     : TentSlab
        outflow      : Optional[Region]
        inflow       : Optional[Region]
        reflect      : Optional[Region]
        """
        ConservationLaw.__init__(self, gfu, tentslab, "burgers", **kwargs)


class Euler(ConservationLaw):
    def __init__(self, gfu, tentslab, **kwargs):
        """
        INPUTS:

        gridfunction : GridFunction
        tentslab     : TentSlab
        outflow      : Optional[Region]
        inflow       : Optional[Region]
        reflect      : Optional[Region]
        """
        ConservationLaw.__init__(self, gfu, tentslab, "euler", **kwargs)


class Wave(ConservationLaw):
    def __init__(self, gfu, tentslab, **kwargs):
        """
        INPUTS:

        gridfunction : GridFunction
        tentslab     : TentSlab
        outflow      : Optional[Region]
        inflow       : Optional[Region]
        reflect      : Optional[Region]
        transparent  : Optional[Region]
        """
        ConservationLaw.__init__(self, gfu, tentslab, "wave", **kwargs)


class Advection(ConservationLaw):
    def __init__(self, gfu, tentslab, **kwargs):
        """
        INPUTS:

        gridfunction : GridFunction
        tentslab     : TentSlab
        outflow      : Optional[Region]
        inflow       : Optional[Region]
        reflect      : Optional[Region]
        transparent  : Optional[Region]
        """
        ConservationLaw.__init__(self, gfu, tentslab, "advection", **kwargs)


class Maxwell(ConservationLaw):
    def __init__(self, gfu, tentslab, **kwargs):
        """
        INPUTS:

        gridfunction : GridFunction
        tentslab     : TentSlab
        outflow      : Optional[Region]
        inflow       : Optional[Region]
        reflect      : Optional[Region]
        """
        ConservationLaw.__init__(self, gfu, tentslab, "maxwell", **kwargs)
