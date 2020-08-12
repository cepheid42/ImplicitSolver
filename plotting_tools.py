import numpy as np
import plotly.graph_objs as go

def get_plot_choice():
    print("1. animated surface.")
    print("2. animated contour.")
    print("3. stacked lines")
    while True:
        choice = int(input('Select plot type: '))

        if 1 <= choice <= 3:
            return choice
        else:
            print('Invalid selection, try again.')


def get_file_choice():
    print("1. ex.")
    print("2. ey.")
    print("3. ez.")
    while True:
        choice = int(input('Select files to plot: '))

        if 1 <= choice <= 3:
            return choice
        else:
            print('Invalid selection, try again.')


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
        filename = f'{folder}/t{q}.csv'
        try:
            temp = np.genfromtxt(filename, delimiter=',', dtype=np.float32).reshape((p.nz, p.ny, p.nx))
            store.append(temp)
            print(f'File {filename}')
        except:
            print(f'Unable to open file {filename}. Continuing.')
            continue

    # Get last file
    try:
        temp = np.genfromtxt(f'{folder}/t{p.nt}.csv', delimiter=',', dtype=np.float32).reshape((p.nz, p.ny, p.nx))
        store.append(temp)
        print(f'File {folder}/t{p.nt}.csv')
    except:
        print(f'Unable to open file {folder}/t{p.nt}.csv. Continuing.')

    return np.asarray(store)


def animated_surface(files, n, p):
    cmin = np.min(files)
    cmax = np.max(files)

    fig = go.Figure()

    for fi in files:
        fig.add_trace(go.Surface(z=fi, opacity=0.7, colorscale='Rainbow', cmin=cmin, cmax=cmax))

    steps = []
    for i in range(len(files)):
        step = dict(method='update', args=[{'visible': [False] * len(files)},
                                           {'title': f'Step {i}'}],)
        step['args'][0]['visible'][i] = True
        steps.append(step)

    sliders = [dict(active=0, currentvalue={'prefix': 'Frequency: '}, pad={'t': 50}, steps=steps)]

    fig.update_layout(scene_aspectmode='cube',
                      scene=dict(xaxis=dict(range=[0, p.nx * p.dx]),
                                 yaxis=dict(range=[0, p.ny * p.dy]),
                                 zaxis=dict(range=[-160, 160])),
                      sliders=sliders)
    fig.show()


def stacked_line(files, n, p):
    lines = []

    for fi in files:
        lines.append(go.Scatter(y=fi))

    fig = go.Figure(data=lines[0])

    for l in range(1, len(lines) - 1):
        fig.add_trace(lines[l])

    fig.show()


def animated_contour(files, n, p):
    zmin = np.min(files)
    zmax = np.max(files)

    fig = go.Figure()

    for fi in files:
        fig.add_trace(go.Contour(z=fi, opacity=0.7, colorscale='Rainbow', zmin=zmin, zmax=zmax))

    steps = []
    for i in range(len(files)):
        step = dict(method='update', args=[{'visible': [False] * len(files)},
                                           {'title': f'Step {i}'}],)
        step['args'][0]['visible'][i] = True
        steps.append(step)

    sliders = [dict(active=0, currentvalue={'prefix': 'Frequency: '}, pad={'t': 50}, steps=steps)]

    fig.update_layout(scene_aspectmode='cube',
                      scene=dict(xaxis=dict(range=[0, p.nx * p.dx]),
                                 yaxis=dict(range=[0, p.ny * p.dy]),
                                 zaxis=dict(range=[-160, 160])),
                      sliders=sliders)
    fig.show()

if __name__ == '__main__':
    params = Params()

    # ex_files = get_files('outputs/ex', params)
    # ey_files = get_files('outputs/ey', params)
    ez_files = get_files('outputs/ez', params)
    # bx_files = get_files('outputs/bx', params)
    # by_files = get_files('outputs/by', params)
    # bz_files = get_files('outputs/bz', params)

    while True:
        plot_choice = get_plot_choice()
        file_choice = get_file_choice()

        if file_choice == 1:
            f = np.asarray([a[params.nz // 2, :, :] for a in ex_files])
            name = 'ex @ nz / 2'
        elif file_choice == 2:
            f = np.asarray([a[params.nz // 2, :, :] for a in ey_files])
            name = 'ey @ nz / 2'
        else:
            f = np.asarray([a[params.nz // 2, :, :] for a in ez_files])
            name = 'ez @ nz / 2'

        if plot_choice == 1:
            animated_surface(f, name, params)
        elif plot_choice == 2:
            animated_contour(f, name, params)
        else:
            stacked_line(f, name, params)

