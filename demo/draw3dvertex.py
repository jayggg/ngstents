from ngstents import TentSlab
from netgen.csg import unit_cube
import ngsolve as ng
import numpy as np
import matplotlib.pyplot as plt

mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=0.3))
resp = input("Enter a mesh vertex number in 0:{}: ".format(mesh.nv))
vnr = int(resp)
tps = TentSlab(mesh, dt=0.5, c=1.0)
ntents = tps.GetNTents()
tents = [(i, tps.GetTent(i)) for i in range(ntents)]
tents = sorted([(i, t) for (i, t) in tents
    if t.vertex == vnr], key=lambda pr: pr[1].level)
plt.ion()
for i,_ in tents:
    plt.close()
    ax = tps.Draw3DTentPlt(i)
    plt.show()
    input("press <Enter> ")



