import numpy as np
import matplotlib.pyplot as plt

# 1D uniform mesh of [0,1]
els = 10*np.array([2**p for p in range(-1, 6)])
h = .1*np.array([2**(-p) for p in range(-1, 6)])
err0 = [0.0, 2.22e-16, 2.2e-16, 1.11e-15, 2.22e-15, 4.44e-15, 8.88e-15]
err1 = [0.0, 6.66e-16, 1.44e-14, 5.77e-14, 2.558e-13, 9.77e-13, 3.83e-12]

plt.loglog(h, err0, 'b.-', label="Geom")
plt.loglog(h, err1, 'r.-', label="H1")
plt.loglog([.1/2**5, .1/2**4], [5.0e-15, 2.5e-15], 'g-', label="$O(h)$")
plt.loglog([.1/2**5, .1/2**4], [2.0e-12, 5e-13], 'c-', label="$O(h^2)$")
plt.title(r"1D unif mesh: increase in grad err with refinement")
plt.xlabel("element size h")
plt.ylabel("max error over all tents")
plt.legend()
plt.show()

# 2D mesh of unit square initially maxh=.2 (obtuse triangles)
h = .2*np.array([2**(-p) for p in range(5)])
err0 = [5.55e-16, 1.22e-15, 1.77e-15, 3.11e-15, 5.77e-15]
err1 = [2.e-15, 2.39e-14, 1.06e-13, 5.19e-13, 2.18e-12]

plt.loglog(h, err0, 'b.-', label="Geom")
plt.loglog(h, err1, 'r.-', label="H1")
plt.loglog([.1/2**3, .1/2**2], [2.5e-15, 1.25e-15], 'g-', label="$O(h)$")
plt.loglog([.1/2**3, .1/2**2], [1.0e-12, 2.5e-13], 'c-', label="$O(h^2)$")
plt.title(r"2D unstruct unif mesh: increase in grad err with refinement")
plt.xlabel("element size h")
plt.ylabel("max error over all tents")
plt.legend()
plt.show()

