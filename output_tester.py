import numpy as np

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
            # print(f'Unable to open file {filename}. Continuing.')
            continue


    return np.asarray(store)

if __name__ == '__main__':
    params = Params()
    first_run = get_files('ez_all_jz_1', params)
    second_run = get_files('outputs/ez', params)

    for i in range(first_run.shape[0]):
        print(np.allclose(first_run[i], second_run[i]))