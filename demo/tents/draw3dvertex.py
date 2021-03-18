"""
Select a mesh vertex by number and draw each tent for which the
specified vertex is the tent vertex
"""

from ngstents import TentSlab
from netgen.csg import unit_cube
import ngsolve as ng
import matplotlib.pyplot as plt

mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=0.3))
resp = input("Enter a mesh vertex number in 0:{}: ".format(mesh.nv-1))
try:
    vnr = int(resp)
except ValueError:
    vnr = -1
if 0 <= vnr < mesh.nv:
    tps = TentSlab(mesh)
    tps.SetMaxWavespeed(1.0)
    tps.PitchTents(dt=0.5)
    ntents = tps.GetNTents()
    tents = [(i, tps.GetTent(i)) for i in range(ntents)]
    tents = sorted([(i, t) for (i, t) in tents
                    if t.vertex == vnr], key=lambda pr: pr[1].level)
    plt.ion()
    for i, _ in tents:
        plt.close()
        ax = tps.Draw3DTentPlt(i)
        plt.show()
        input("press <Enter> ")
else:
    print("invalid vertex")
