"""
Select a 3D tents to visualize in pyplot
"""

from ngstents import TentSlab
from netgen.csg import unit_cube
import ngsolve as ng
import matplotlib.pyplot as plt

mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=0.3))
print("creating tent pitched slab")
tps = TentSlab(mesh)
tps.SetMaxWavespeed(1.0)
tps.PitchTents(dt=0.5)
ntents = tps.GetNTents()
print("The slab has {} tents.".format(ntents))
resp = ""
plt.ion()
while True:
    resp = input("Enter a tent number in 0:{} or q to quit: ".format(ntents-1))
    plt.close()
    if resp.lower() == 'q':
        break
    try:
        n = int(resp)
    except ValueError:
        continue
    if 0 <= n < ntents:
        ax = tps.Draw3DTentPlt(n)
        plt.show()
    else:
        print("Invalid tent number")
