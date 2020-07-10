import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Params:
    def __init__(self):
        self.nx = 0
        self.ny = 0
        self.nz = 0

        self.dx = 0
        self.dy = 0
        self.dz = 0

        self.nt = 0
        self.step = 0
        self.load_params()

    def load_params(self):
        with open('outputs/params.csv', 'r') as p:
            lines = p.readlines()
            self.nx, self.ny, self.nz = [int(val.strip()) for val in lines[0].split(',')]
            self.dx, self.dy, self.dz = [float(val.strip()) for val in lines[1].split(',')]
            self.nt = int(lines[2].strip())
            self.step = int(lines[3].strip())


def get_files(folder, params):
    store = []

    for q in range(0, params.nt, params.step):
        filename = f'outputs/{folder}/t{q}.csv'
        try:
            temp = np.genfromtxt(filename, dtype=np.float64).reshape((params.nz, params.ny, params.nx))
            store.append(temp)
        except:
            print(f'Unable to open file {filename}. Continuing.')

    return np.asarray(store)


if __name__ == '__main__':
    p = Params()

    fig = plt.figure(figsize=(20, 16))
    ax1 = fig.add_subplot(131, projection='3d')
    ax1.set_title('Ex')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Time (line)')

    ax2 = fig.add_subplot(132, projection='3d')
    ax2.set_title('Ey')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Time (line)')

    ax3 = fig.add_subplot(133, projection='3d')
    ax3.set_title('Ex')
    ax3.set_xlabel('X')
    ax3.set_ylabel('Time (line)')

    ex_files = get_files('ex', p.nt, p.step)
    ey_files = get_files('ey', p.nt, p.step)

    x_range = np.arange(p.nx - 1)

    count1 = 0
    for f in ex_files:
        q = count1 * p.step
        ax1.plot(x_range, f[1:] + q, 'k', lw=1, zorder=p.nt - q)
        ax1.fill_between(x_range, f[1:] + q, facecolor='w', lw=0, zorder=p.nt - q - 1)
        count1 += 1

    count2 = 0
    for f in ey_files:
        q = count1 * p.step
        ax2.plot(x_range, f[1:] + q, 'k', lw=1, zorder=p.nt - q)
        ax2.fill_between(x_range, f[1:] + q, facecolor='w', lw=0, zorder=p.nt - q - 1)
        count2 += 1

    count3 = 0
    for f in ez_files:
        q = count2 * p.step
        ax3.plot(x_range, f[1:] + q, 'k', lw=1, zorder=p.nt - q)
        ax3.fill_between(x_range, f[1:] + q, facecolor='w', lw=0, zorder=p.nt - q - 1)
        count3 += 1

    plt.show()
