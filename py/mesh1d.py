from netgen.meshing import Mesh, MeshPoint, Element0D, Element1D
from netgen.csg import Pnt
import matplotlib.pyplot as plt

def Make1DMeshSimple(nels, save=False, filename="mesh1d.vol"):
    """
    Create a uniform 1D mesh of the interval [0, 1] with inflow and
    outflow boundary conditions with the specified number of elements
    """
    m = Mesh()
    m.dim = 1

    pnums = []
    for i in range(0, nels+1):
        pnums.append(m.Add(MeshPoint(Pnt(i/nels, 0, 0))))

    for i in range(0, nels):
        m.Add(Element1D([pnums[i], pnums[i+1]], index=1))

    m.Add(Element0D(pnums[0], index=1))
    m.Add(Element0D(pnums[nel], index=3))
    if(save):
        m.Save(filename)
    else:
        return m


def Make1DMeshSpecified(pts, bc=[1, 3], bcname=None,
                        save=False, filename="mesh1d.vol"):
    """
    Create a 1D mesh of the interval [0,1] with inflow and outflow
    boundary conditions with specified verices (not necessarily uniformly
    spaced)
    """
    m = Mesh()
    m.dim = 1
    pnums = []
    for pt in pts:
        pnums.append(m.Add(MeshPoint(Pnt(pt, 0, 0))))
    for j in range(len(pts)-1):
        m.Add(Element1D([pnums[j], pnums[j+1]], index=1))
    if(bcname):
        m.Add(Element0D(pnums[0], index=1))
        m.Add(Element0D(pnums[-1], index=2))
        if(isinstance(bcname, list)):
            m.SetBCName(0, bcname[0])
            m.SetBCName(1, bcname[1])
        else:
            m.SetBCName(0, bcname)
            m.SetBCName(1, bcname)
    else:
        m.Add(Element0D(pnums[0], index=bc[0]))
        m.Add(Element0D(pnums[-1], index=bc[1]))

    if(save):
        m.Save(filename)
    else:
        return m


def Make1DMesh(domains, nels, bc=[1, 3], bcname=None,
               save=False, filename="mesh1d.vol"):
    """
    Generate a 1D mesh consisting of multiple domains with specified endpoints
    and boundary conditions.  Named boundary conditions may optionally be set.

    domains: nested list containing left and right endpoints for each domain
    """ 
    m = Mesh()
    m.dim = 1
    pnums = []
    start = 0
    for i in range(0, len(domains)):
        xmin, xmax = domains[i][0], domains[i][1]
        nel = nels[i]
        npnts = len(pnums)-start
        for j in range(start, nel+1):
            pnums.append(m.Add(MeshPoint(Pnt(xmin+(xmax-xmin)*j/nel, 0, 0))))

        for j in range(npnts, npnts+nel):
            m.Add(Element1D([pnums[j], pnums[j+1]], index=i+1))
        if(i == 0):
            start = 1

    if(bcname):
        m.Add(Element0D(pnums[0], index=bc[0]))
        m.Add(Element0D(pnums[-1], index=bc[1]))
        if(isinstance(bcname, list)):
            m.SetBCName(0, bcname[0])
            m.SetBCName(1, bcname[1])
        else:
            m.SetBCName(0, bcname)
            m.SetBCName(1, bcname)
    else:
        m.Add(Element0D(pnums[0], index=bc[0]))
        m.Add(Element0D(pnums[-1], index=bc[1]))

    if(save):
        m.Save(filename)
    else:
        return m


def Make1DPeriodicMesh(nel):
    """
    Generate a 1D periodic mesh of the interval [0,1] with inflow and outflow
    boundary conditions
    """
    m = Mesh()
    m.dim = 1
    pnums = []
    for i in range(0, nel+1):
        pnums.append(m.Add(MeshPoint(Pnt(i/nel, 0, 0))))

    for i in range(0, nel):
        m.Add(Element1D([pnums[i], pnums[i+1]], index=1))

    m.Add(Element0D(pnums[0], index=1))
    m.Add(Element0D(pnums[nel], index=3))
    m.AddPointIdentification(pnums[0], pnums[nel], identnr=1, type=2)
    return m
