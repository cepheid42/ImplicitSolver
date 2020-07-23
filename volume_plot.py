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


def get_files(folder, p):
    store = []

    for q in range(0, p.nt, p.step):
        filename = f'{folder}/t{q}.csv'
        try:
            temp = np.genfromtxt(filename, delimiter=',', dtype=np.float32).reshape((p.nz, p.ny, p.nx))
            store.append(temp)
        except:
            print(f'Unable to open file {filename}. Continuing.')
            continue


    return np.asarray(store)

def animated_surface(files, p):
    frames = [go.Frame(data=go.Surface(z=files[i][:, p.ny // 2, :])) for i in range(1, len(files) - 1)]

    fig = go.Figure(
        data=go.Surface(z=files[0][:, p.ny // 2, :]),
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

def stacked_line(files, p):
    lines = []

    for f in files:
        lines.append(go.Scatter(y=f[p.nz // 2, p.ny // 2, :]))

    fig = go.Figure(data=lines[0])

    for l in range(1, len(lines) - 1):
        fig.add_trace(lines[l])

    fig.show()


if __name__ == '__main__':
    params = Params()
    new_ez = 'outputs/ez'
    old_ez1 = 'ez_all_jz_1'
    old_ez2 = 'ez_all_jz_1_then_0'

    ez_files = get_files(new_ez, params)

    stacked_line(ez_files, params)
