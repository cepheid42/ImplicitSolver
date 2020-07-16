import numpy as np
import matplotlib.pyplot as plt

nx = 97
nt = 138
step = 10


def get_files(folder):
    store = []

    for qq in range(0, nt, step):
        filename = f'outputs/{folder}/t{qq}.csv'
        try:
            temp = np.genfromtxt(filename, delimiter=',', dtype=np.float32)
            store.append(temp)
        except:
            print(f'Unable to open file {filename}. Continuing.')

    return np.asarray(store)


if __name__ == '__main__':
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ez_files = get_files('ez')

    x_range = np.arange(nx - 1)
    # x = np.linspace(0, 100, 1000)
    # y = np.sin(x)
    #
    # ax.plot(x, y)

    count = 0
    for f in ez_files:
        q = count * step
        ax.plot(x_range, f[1:] + q, 'k', lw=1, zorder=nt - q)
        ax.fill_between(x_range, f[1:] + q, facecolor='w', lw=0, zorder=nt - q - 1)
        count += 1

    plt.show()
