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
            temp = np.genfromtxt(filename, delimiter=',', dtype=np.float32).reshape((params.nx, params.ny, params.nz))
            store.append(temp)
        except:
            print(f'Unable to open file {filename}. Continuing.')


    return np.asarray(store)


if __name__ == '__main__':
    p = Params()
    ey_files = get_files('ey', p)

    frames = [go.Frame(data=go.Surface(z=ey_files[i][p.nx // 2, :, :])) for i in range(1, len(ey_files))]


    fig = go.Figure(
        data=go.Surface(z=ey_files[0][p.nx // 2, :, :]),
        layout=go.Layout(
            xaxis=dict(autorange=True),
            yaxis=dict(autorange=True),
            title="Ey",
            updatemenus=[dict(
                type="buttons",
                buttons=[dict(label="Play",
                              method="animate",
                              args=[None])])]
        ),
        frames=frames
    )



    fig.show()
