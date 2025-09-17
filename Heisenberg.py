import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

filein = open("tmax.txt")
tmax = float(filein.readline())
Ntraj = int(filein.readline())

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

rates = read_file_columns("rates.txt")
trace = read_file_columns("trace.txt")
obs = read_file_columns("obs.txt")
obs_avg = read_file_columns("obs_avg.txt")
traj_x = read_file_columns("traj_x.txt")
traj_z = read_file_columns("traj_z.txt")

t = np.linspace(0, tmax, len(obs[0]))
markevery = int(len(obs[0])/30)

fig, ax = plt.subplots(1,1, figsize=(9,5), sharex=True, sharey=False)
fontsize = 16

ax.plot(t, obs[0], '--', color="red", label=r'$x$', marker='o', fillstyle='none', markevery=len(obs[0]))
ax.plot(t, obs[1], color="blue", label=r'$z$', marker='x', markevery=len(obs[1]))
ax.plot(t, obs_avg[0], color="red", marker='o', fillstyle='none', markevery=markevery, linewidth=0)
ax.plot(t, obs_avg[1], color="blue", marker='x', fillstyle='none', markevery=markevery, linewidth=0)
ax.legend(loc="lower right", fontsize=fontsize-2)
ax.set_ylabel(r'$\operatorname{tr}[X(t)\sigma_\alpha]$', fontsize=fontsize)
ax.set_xlabel(r'$t$', fontsize=fontsize)
for i in range(Ntraj):
    tx = np.array(traj_x[i])
    tz = np.array(traj_z[i])
    mask_x = tx > -10
    mask_z = tz > -10
    ax.plot(t[mask_x], tx[mask_x], '--', color="red", alpha=0.1)
    ax.plot(t[mask_z], tz[mask_z], color="blue", alpha=0.1)

size_inset_rates = .3
size_inset_tr = .3

axx = ax.inset_axes([.05,.02,size_inset_rates,size_inset_rates])
axx.plot(t, rates[0], label=r'$\gamma_-$')
axx.plot(t, rates[1], label=r'$\gamma_+$')
axx.plot(t, rates[2], label=r'$\epsilon$')
#axx.axhline(0, color='black', linewidth=1.)
axx.legend(loc = "upper right")
axx.set_xticks([])
axx.set_yscale('log')

axx = ax.inset_axes([.98-size_inset_tr,.98-size_inset_tr,size_inset_tr,size_inset_tr])
axx.plot(t, trace[0], color='green')
axx.plot(t, trace[1], color='green', marker='d', fillstyle='none', markevery=markevery*3, linewidth=0)
axx.set_ylabel(r'$\operatorname{tr}[X(t)]$', fontsize=fontsize-3)
axx.set_xticks([])

plt.tick_params(axis='both', which='major', labelsize=fontsize-3)
plt.savefig("Heis_TNP.png", dpi=300)
plt.show()