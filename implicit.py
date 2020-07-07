import numpy as np
import matplotlib.pyplot as plt
from time import time
from numba import jit

# Constants
c0 = 299792458.0  # speed of light, m/s
eps0 = 8.85418782e-12
mu0 = 1.25663706e-6

freq = 14.0e6
lamb = c0 / freq
T = 1.0 / freq
t0 = 6.0 * T
sig0 = 2.0 * T
n0 = 2.0

x_resolution = 64
num_wavelengths = 20

nx = num_wavelengths * x_resolution + 1
dx = (num_wavelengths * lamb - 0.0) / (nx - 1)

# courant number
C = 1.0

dt = C * dx / c0
nt = int(num_wavelengths * lamb / c0 / dt + 1)

# [array, col, row]
Ez = np.zeros(nx)

By = np.zeros(nx)

J = np.zeros(nx)

ez_nhalf = np.zeros(nx)
ez_nhalf_rhs = np.zeros(nx)
ez_one = np.zeros(nx)
ez_one_rhs = np.zeros(nx)

coef1 = 1.0 / (8 * eps0 * mu0) * (dt / dx) ** 2
a_0 = - coef1 * np.ones(nx - 1)
b_0 = 0.5 - coef1 * np.full(nx, -2.0)
c_0 = - coef1 * np.ones(nx - 1)


@jit(nopython=True)
def TDMAsolver(a, b, c, d):
    n = len(d)  # number of equations
    ac = np.copy(a)
    bc = np.copy(b)
    cc = np.copy(c)
    dc = np.copy(d)

    for it in range(1, n):
        mc = ac[it - 1] / bc[it - 1]
        bc[it] = bc[it] - mc * cc[it - 1]
        dc[it] = dc[it] - mc * dc[it - 1]

    xc = bc
    xc[-1] = dc[-1] / bc[-1]

    for il in range(n - 2, -1, -1):
        xc[il] = (dc[il] - cc[il] * xc[il + 1]) / bc[il]

    return xc


print(f'Running {nt} iterations.')
start = time()

fig = plt.figure()
ax1 = fig.add_subplot(111)

for q in range(nt):
    t = q * dt
    # Source
    J[0] = 9.89399 * np.sin(2.0 * np.pi * freq * (t + 0.5 * dt)) * np.exp(-1 * ((t + 0.5 * dt) - t0) ** n0 / (2.0 * sig0 ** n0))

    """ N -> N + 1/2 """
    # Implicit update to auxillary
    ez_nhalf_rhs[1:] = Ez[1:] + (1 / eps0) * (0.5 * dt / dx) * (By[1:] - By[:nx - 1]) - (0.5 * dt / eps0) * J[:nx - 1]

    ez_nhalf = TDMAsolver(a_0, b_0, c_0, ez_nhalf_rhs)

    Ez[:nx - 1] = ez_nhalf[:nx - 1] - Ez[:nx - 1]

    # Explicit update magnetic field
    By[:nx - 1] = By[:nx - 1] + (1 / mu0) * (0.5 * dt / dx) * (ez_nhalf[1:] - ez_nhalf[:nx - 1])

    """ N + 1/2 -> N + 1"""
    # Implicit update to auxillary
    ez_one_rhs = Ez
    ez_one = 2.0 * ez_one_rhs

    Ez = ez_one - Ez

    if q % 50 == 0:
        print(f'{q}/{nt}: {time() - start:.3f}s')
        ax1.plot(np.arange(nx - 1), Ez[1:] + q, 'k', lw=1, zorder=nt - q)
        ax1.fill_between(np.arange(nx - 1), Ez[1:] + q, facecolor='w', lw=0, zorder=nt - q - 1)

plt.show()
