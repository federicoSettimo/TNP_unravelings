import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

filein = open("tmax.txt")
tmax = float(filein.readline())
z0 = float(filein.readline())
z1 = float(filein.readline())
z2 = float(filein.readline())
z3 = float(filein.readline())
z4 = float(filein.readline())
filein.close()

def read_file_columns(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            numbers = line.split()
            for i, number in enumerate(numbers):
                if len(data) <= i:
                    data.append([])
                data[i].append(float(number))
    return np.array(data)

mom = read_file_columns("moments.txt")
mom_ex = read_file_columns("moments_ex.txt")

trace_0 = read_file_columns("trace_0.txt")
trace_ex_0 = read_file_columns("trace_ex_0.txt")
trace_1 = read_file_columns("trace_1.txt")
trace_ex_1 = read_file_columns("trace_ex_1.txt")
trace_2 = read_file_columns("trace_2.txt")
trace_ex_2 = read_file_columns("trace_ex_2.txt")
trace_3 = read_file_columns("trace_3.txt")
trace_ex_3 = read_file_columns("trace_ex_3.txt")
trace_4 = read_file_columns("trace_4.txt")
trace_ex_4 = read_file_columns("trace_ex_4.txt")

t = np.linspace(0, tmax, len(mom[0]))

fig, ax = plt.subplots(1,1, figsize=(9,5))
fontsize = 16
size_inset = .35
axx = ax.inset_axes([.04,.95 - size_inset,size_inset,size_inset])
axx.set_yticks([])
axx = axx.twinx()

markevery = len(t)//20
ax.plot(t[1:-1], mom[0][1:-1], color='red', markevery=markevery, marker='o', lw=0, fillstyle='none')
ax.plot(t[1:-1], mom_ex[0][1:-1], color='red', label=f'$\mu_1$', markevery=len(t), fillstyle='none', marker='o')
ax.plot(t[1:-1], mom[1][1:-1], color='green', markevery=markevery, marker='x', lw=0, fillstyle='none')
ax.plot(t[1:-1], mom_ex[1][1:-1], color='green', label=f'$\mu_2$', markevery=len(t), fillstyle='none', marker='x')
ax.plot(t[1:-1], mom[2][1:-1], color='blue', markevery=markevery, marker='d', lw=0, fillstyle='none')
ax.plot(t[1:-1], mom_ex[2][1:-1], color='blue', label=f'$\mu_3$', markevery=len(t), fillstyle='none', marker='d')
ax.plot(t[1:-1], mom[3][1:-1], color='orange', markevery=markevery, marker='s', lw=0, fillstyle='none')
ax.plot(t[1:-1], mom_ex[3][1:-1], color='orange', label=f'$\mu_4$', markevery=len(t), fillstyle='none', marker='s')


ax.legend(loc="lower right", fontsize=fontsize)
#ax.set_ylabel(r'Moments', fontsize=fontsize)
ax.set_xlabel(r'$t$', fontsize=fontsize)
#ax.set_ylim([-.55,1.55])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
markevery = len(t)//10
axx.plot(t, trace_0[0], color=colors[0], markevery=markevery, marker='o', lw=0, fillstyle='none')
axx.plot(t, trace_ex_0[0], label=rf'$\zeta={z0}$', color=colors[0], markevery=len(t), fillstyle='none', marker='o')
axx.plot(t, trace_1[0], color=colors[1], markevery=markevery, marker='x', lw=0, fillstyle='none')
axx.plot(t, trace_ex_1[0], label=rf'$\zeta={z1}$', color=colors[1], markevery=len(t), fillstyle='none', marker='x')
axx.plot(t, trace_2[0], color=colors[2], markevery=markevery, marker='d', lw=0, fillstyle='none')
axx.plot(t, trace_ex_2[0], label=rf'$\zeta={z2}$', color=colors[2], markevery=len(t), fillstyle='none', marker='d')
axx.plot(t, trace_3[0], color=colors[3], markevery=markevery, marker='s', lw=0, fillstyle='none')
axx.plot(t, trace_ex_3[0], label=rf'$\zeta={z3}$', color=colors[3], markevery=len(t), fillstyle='none', marker='s')
axx.plot(t, trace_4[0], color=colors[4], markevery=markevery, marker='v', lw=0, fillstyle='none')
axx.plot(t, trace_ex_4[0], label=rf'$\zeta={z4}$', color=colors[4], markevery=len(t), fillstyle='none', marker='v')

axx.set_ylabel(r'$1-\operatorname{tr}\rho_\zeta$', fontsize=fontsize - 2)
axx.set_xticks([])


plt.savefig('photon_counting_moments.png', dpi=300)
plt.show()