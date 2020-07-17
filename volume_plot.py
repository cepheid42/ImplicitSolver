import numpy as np
import plotly.graph_objs as go
import plotly.express as px

class Params:
    def __init__(self):
        self.nx = 0
        self.ny = 0
        self.nz = 0

        self.dx = 0
        self.dy = 0
        self.dz = 0

        self.nt = 0
        self.dt = 0
        self.step = 0

        self.load_params()

    def load_params(self):
        with open('outputs/params.csv', 'r') as p:
            lines = p.readlines()
            self.nx, self.ny, self.nz = [int(val.strip()) for val in lines[0].split(',')]
            self.dx, self.dy, self.dz = [float(val.strip()) for val in lines[1].split(',')]
            self.nt = int(lines[2].split(',')[0].strip())
            self.dt = float(lines[2].split(',')[1].strip())
            self.step = int(lines[3].strip())


def get_files(folder, params):
    store = []

    for q in range(0, params.nt, params.step):
        filename = f'outputs/{folder}/t{q}.csv'
        try:
            temp = np.genfromtxt(filename, delimiter=',', dtype=np.float32).reshape((params.nz, params.ny, params.nx))
            store.append(temp)
        except:
            print(f'Unable to open file {filename}. Continuing.')


    return np.asarray(store)

def animated_surface():
    frames = [go.Frame(data=go.Surface(z=ez_files[i][1:p.nx, p.ny // 2, 1:p.nz])) for i in range(1, len(ez_files) - 1)]

    fig = go.Figure(
        data=go.Surface(z=ez_files[0][1:p.nx, p.ny // 2, 1:p.nz]),
        layout=go.Layout(
            xaxis=dict(autorange=True),
            yaxis=dict(autorange=True),
            title="ez",
            updatemenus=[dict(
                type="buttons",
                buttons=[dict(label="Play",
                              method="animate",
                              args=[None])])]
        ),
        frames=frames
    )

    fig.show()


if __name__ == '__main__':
    p = Params()
    ez_files = get_files('ez', p)

    x = np.arange(p.nt * p.dt)

    y = [ez_files[i][5, 15, 25] for i in range(len(ez_files))]

    fig = go.Figure(data=go.Scatter(x=x, y=y))

    fig.show()
