import numpy as np
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import matplotlib.animation as animation


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
    print(f'Reading from "{folder}"')
    store = []

    for q in range(0, p.nt, p.step):
        print(f'File "t{q}.csv"')
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

def surface_plot(data, p):
    fig = go.Figure(data=go.Surface(z=data))
    fig.show()

def stacked_line(name, files, p):
    lines = []

    for f in files:
        lines.append(go.Scatter(y=f[p.nz // 2, p.ny // 2, :]))

    fig = go.Figure(data=lines[0])

    for l in range(1, len(lines) - 1):
        fig.add_trace(lines[l])

    fig.update_layout(title=name)
    fig.show()

def line_plot_active(folder, p):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    def animate(i):
        filename = f'{folder}/t{i}.csv'
        yar = []
        while True:
            try:
                dataArray = np.genfromtxt(filename, delimiter=',', dtype=np.float32).reshape((p.nz, p.ny, p.nx))
                yar = dataArray[p.nz // 2, p.ny // 2, :]
            except:
                continue
            break
        xar = np.arange(p.nx)
        ax1.clear()
        ax1.plot(xar, yar)
        ax1.annotate(f'Time Step: {i}', xy=(0.7, 0.9), xycoords='axes fraction', annotation_clip=False)

    anim = animation.FuncAnimation(fig, animate, interval=1000)
    plt.show()

if __name__ == '__main__':
    params = Params()

    ex_files = get_files('outputs/ex', params)
    ey_files = get_files('outputs/ey', params)
    ez_files = get_files('outputs/ez', params)
    bx_files = get_files('outputs/bx', params)
    by_files = get_files('outputs/by', params)
    bz_files = get_files('outputs/bz', params)

    stacked_line("ex", ex_files, params)
    stacked_line("ey", ey_files, params)
    stacked_line("ez", ez_files, params)
    stacked_line("bx", bx_files, params)
    stacked_line("by", by_files, params)
    stacked_line("bz", bz_files, params)

    # surface_plot(ez_files[0][:, params.ny // 2, :], params)
    # animated_surface(ez_files, params)
    # line_plot_active(new_ez, params)
