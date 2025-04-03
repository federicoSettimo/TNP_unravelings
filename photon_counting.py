import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

filein = open("tmax.txt")
tmax = float(filein.readline())

def read_file_columns(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            numbers = line.split()
            for i, number in enumerate(numbers):
                if len(data) <= i:
                    data.append([])
                data[i].append(float(number))
    return data

trace_0 = read_file_columns("trace_0.txt")
obs_0 = read_file_columns("obs_0.txt")
obs_avg_0 = read_file_columns("obs_avg_0.txt")
trace_p = read_file_columns("trace_+.txt")
obs_p = read_file_columns("obs_+.txt")
obs_avg_p = read_file_columns("obs_avg_+.txt")
trace_m = read_file_columns("trace_-.txt")
obs_m = read_file_columns("obs_-.txt")
obs_avg_m = read_file_columns("obs_avg_-.txt")

t = np.linspace(0, tmax, len(obs_0[0]))
markevery = int(len(obs_0[0])/30)

fig, ax = plt.subplots(1,1, figsize=(9,5), sharex=True, sharey=False)

if False:
    ax.plot(t, obs_0[0], color="red", label=r'$\zeta=0$', marker='o', fillstyle='none', markevery=len(obs_0[0]))
    ax.plot(t, obs_0[1], '--', color="red", marker='o', fillstyle='none', markevery=len(obs_0[0]))
    ax.plot(t, obs_avg_0[0], color="red", marker='o', fillstyle='none', markevery=markevery, linewidth=0)
    ax.plot(t, obs_avg_0[1], color="red", marker='o', fillstyle='none', markevery=markevery, linewidth=0)

    ax.plot(t, obs_p[0], color="green", label=r'$\zeta=.01$', marker='x', fillstyle='none', markevery=len(obs_0[0]))
    ax.plot(t, obs_p[1], '--', color="green", marker='x', fillstyle='none', markevery=len(obs_0[0]))
    ax.plot(t, obs_avg_p[0], color="green", marker='x', fillstyle='none', markevery=markevery, linewidth=0)
    ax.plot(t, obs_avg_p[1], color="green", marker='x', fillstyle='none', markevery=markevery, linewidth=0)

    ax.plot(t, obs_m[0], color="blue", label=r'$\zeta=-1$', marker='d', fillstyle='none', markevery=len(obs_0[0]))
    ax.plot(t, obs_m[1], '--', color="blue", marker='d', fillstyle='none', markevery=len(obs_0[0]))
    ax.plot(t, obs_avg_m[0], color="blue", marker='d', fillstyle='none', markevery=markevery, linewidth=0)
    ax.plot(t, obs_avg_m[1], color="blue", marker='d', fillstyle='none', markevery=markevery, linewidth=0)

    ax.legend(loc="upper left")
    ax.set_ylabel(r'$\langle i\vert\rho_\zeta\vert i\rangle$')
    ax.set_xlabel(r'$t$')
    ax.set_title('Photon counting master equation')

    size_inset = .2

    axx = ax.inset_axes([.98-size_inset,.91-size_inset,size_inset,size_inset])
    axx.plot(t, trace_0[0], color='red')
    axx.plot(t, trace_0[1], color='red', marker='o', fillstyle='none', markevery=markevery*3, linewidth=0)
    axx.plot(t, trace_p[0], color='green')
    axx.plot(t, trace_p[1], color='green', marker='x', fillstyle='none', markevery=markevery*3, linewidth=0)
    axx.plot(t, trace_m[0], color='blue')
    axx.plot(t, trace_m[1], color='blue', marker='d', fillstyle='none', markevery=markevery*3, linewidth=0)
    axx.set_ylabel(r'$\operatorname{tr}[\rho_\zeta]$')
    axx.set_xticks([])
else:
    ax.plot(t, trace_0[0], color='red', label=r'$\zeta=0$', marker='o', fillstyle='none', markevery=len(obs_0[0]))
    ax.plot(t, trace_0[1], color='red', marker='o', fillstyle='none', markevery=markevery, linewidth=0)
    ax.plot(t, trace_p[0], color='green', label=r'$\zeta=.01$', marker='x', fillstyle='none', markevery=len(obs_0[0]))
    ax.plot(t, trace_p[1], color='green', marker='x', fillstyle='none', markevery=markevery, linewidth=0)
    ax.plot(t, trace_m[0], color='blue', label=r'$\zeta=-1$', marker='d', fillstyle='none', markevery=len(obs_0[0]))
    ax.plot(t, trace_m[1], color='blue', marker='d', fillstyle='none', markevery=markevery, linewidth=0)
    ax.set_ylabel(r'$\operatorname{tr}[\rho_\zeta]$')
    ax.legend(loc="lower left")
    ax.set_xlabel(r'$t$')
    ax.set_title('Photon counting master equation')

plt.show()