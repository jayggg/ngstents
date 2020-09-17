import numpy as np
import numpy.linalg as la
from pytest import approx

from erickson import (calc_pH, calc_gradHtau, closest_pt_to_Trig)

def test_calc_pH():
    pts = np.array([[1., 0, 0], [0, 1, 0], [0, 0, 1]])
    po = np.array([2, 1, 3.])
    pH = calc_pH(pts, po)
    assert pH[2] + pH[0] + pH[1] == approx(1)


def test_calcgradHtau():
    pts = np.array([[1, 0, 0, 1.0], [0, 1, 0, 0], [0, 0, 1, 0]])
    mu, nu, v, n = calc_gradHtau(pts)
    print(mu*v + nu*n)
    assert mu*v + nu*n == approx(np.array([2/3., -1/3, -1/3]))
    pts = np.array([[1, 0, 0, 1.0], [0, 1, 0, 1.0], [0, 0, 1, 1.0]])
    mu, nu, v, n = calc_gradHtau(pts)
    print(mu*v + nu*n)
    assert mu*v + nu*n == approx(np.zeros(3))
    assert np.sqrt(mu*mu + nu*nu) == approx(0.)


def test_closestpt():
    trig = np.array([[0., 0, 0], [3, 0, 0], [1, 1, 0]])
    pts = np.array([[-1., -1, 0], [0, 4/3, 0], [1, 4/3, 0],
                    [2, 1, 0], [4, 0, 0], [1, .5, 0]])
    cps = np.array([[0., 0, 0], [2/3, 2/3, 0], [1, 1, 0],
                    [1.8, .6, 0], [3, 0, 0], [1, .5, 0]])
    for p, cp in zip(pts, cps):
        assert closest_pt_to_Trig(trig, p) == approx(cp)
