import numpy as np
import plotly.graph_objs as go

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
    ez_files = get_files('ez', p)
    ez = ez_files[0]

    ix = p.nx * p.dx
    iy = p.ny * p.dy
    iz = p.nz * p.dz

    X, Y, Z = np.mgrid[0:iz:p.nz, 0:iy:p.ny, 0:ix:p.nx]

    fig = go.Figure(data=go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=ez.flatten(),
        isomin=0.1,
        isomax=0.8,
        opacity=0.1, # needs to be small to see through all surfaces
        surface_count=17, # needs to be a large number for good volume rendering
    ))
    fig.show()
