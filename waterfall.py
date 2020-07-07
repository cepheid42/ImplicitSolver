import numpy as np
import matplotlib.pyplot as plt

step = 10
nt = 1281
nx = 1281

def get_files():
    store = []

    for f in range(0, nt, step):
        filename = f'outputs/ez_q{f}.csv'
        try:
            temp = np.genfromtxt(filename, dtype=np.float64)
            store.append(temp)
        except:
            print(f'Unable to open file {filename}. Continuing.')

    return np.asarray(store)

if __name__ == '__main__':
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    x_range = np.arange(nx - 1)

    files = get_files()
    count = 1
    ax1.plot(x_range, files[35][1:], 'k', lw=1)
    # for f in files:
    #     q = count * step
    #     ax1.plot(x_range, f[1:] + q, 'k', lw=1, zorder=nt - q)
    #     ax1.fill_between(x_range, f[1:] + q, facecolor='w', lw=0, zorder=nt - q - 1)
    #     count += 1

    plt.show()